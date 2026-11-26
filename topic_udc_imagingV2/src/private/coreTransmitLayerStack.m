function transmissionAmplitudeByRegion = coreTransmitLayerStack(wavelengths, layerStack, regionPatterns, layerImageDir)
%coreTransmitLayerStack 求各区域法向复透射，以入射功率归一化振幅。
    validateattributes(layerStack.incidentN, {'numeric'}, {'scalar','real','finite','positive'});
    validateattributes(layerStack.exitN, {'numeric'}, {'scalar','real','finite','positive'});
    layerNames = string(layerStack.layerNames);
    backgroundMode = "opaque";
    if isfield(layerStack,'backgroundMode'), backgroundMode = string(layerStack.backgroundMode); end
    mustBeMember(backgroundMode,["opaque","medium"]);
    usesSpecifiedTransmittance = false(1,numel(layerNames));
    specifiedAmplitude = zeros(numel(layerNames),numel(wavelengths));
    complexRefractiveIndex = zeros(numel(layerNames),numel(wavelengths));
    layerThickness = zeros(1,numel(layerNames));
    for j = 1:numel(layerNames)
        layer = layerStack.layers.(layerNames(j));
        mustBeMember(string(layer.method), ["transmittance","nk"]);
        usesSpecifiedTransmittance(j) = string(layer.method) == "transmittance";
        if usesSpecifiedTransmittance(j)
            validateattributes(layer.transmittance, {'numeric'}, {'vector','real','finite','>=',0,'<=',1});
            assert(isscalar(layer.transmittance) || numel(layer.transmittance)==numel(wavelengths), ...
                'UDCImagingV2:TransmittanceSize', '透过率应为标量或与波长一一对应。');
            specifiedAmplitude(j,:) = sqrt(reshape(layer.transmittance,1,[]));
        else
            validateattributes(layer.thickness, {'numeric'}, {'scalar','real','finite','nonnegative'});
            complexRefractiveIndex(j,:) = layerIndex(layer,wavelengths,layerImageDir);
            layerThickness(j) = layer.thickness;
        end
    end
    referenceHeight = 0;
    for level = unique(layerStack.layerLevels)
        referenceHeight = referenceHeight + max(layerThickness(layerStack.layerLevels==level));
    end

    transmissionAmplitudeByRegion = complex(zeros(size(regionPatterns,1),numel(wavelengths)));
    for region = 1:size(regionPatterns,1)
        presentLayers = find(regionPatterns(region,:));
        if isempty(presentLayers) && backgroundMode=="opaque"
            continue
        end
        specifiedLayers = presentLayers(usesSpecifiedTransmittance(presentLayers));
        if ~isempty(specifiedLayers)
            assert(isscalar(presentLayers), 'UDCImagingV2:DirectOverlap', ...
                '给定功率透过率适用于整个材料区域，该区域不可与其他透射层重叠。');
            transmissionAmplitudeByRegion(region,:) = specifiedAmplitude(specifiedLayers,:);
            continue
        end
        % 图层从下到上存储，从出射界面向入射侧递推。
        filmLayers = find(regionPatterns(region,:) & layerThickness>0);
        mediumRefractiveIndices = [layerStack.exitN*ones(size(wavelengths)); ...
            complexRefractiveIndex(filmLayers,:); layerStack.incidentN*ones(size(wavelengths))];
        reflectionAmplitude = zeros(size(wavelengths));
        transmissionAmplitude = ones(size(wavelengths));
        filmPropagationFactor = ones(size(wavelengths));
        for j = 2:size(mediumRefractiveIndices,1)
            lowerIndex = mediumRefractiveIndices(j-1,:);
            upperIndex = mediumRefractiveIndices(j,:);
            interfaceReflection = (upperIndex-lowerIndex)./(upperIndex+lowerIndex);
            reflectionDenominator = 1 + interfaceReflection.*reflectionAmplitude.*filmPropagationFactor.^2;
            transmissionAmplitude = (1+interfaceReflection).*filmPropagationFactor.*transmissionAmplitude ...
                ./reflectionDenominator;
            reflectionAmplitude = (interfaceReflection+reflectionAmplitude.*filmPropagationFactor.^2) ...
                ./reflectionDenominator;
            if j <= numel(filmLayers)+1
                filmPropagationFactor = exp(2i*pi*upperIndex*layerThickness(filmLayers(j-1))./wavelengths);
            end
        end
        referencePlaneGap = max(0,referenceHeight-sum(layerThickness(filmLayers)));
        transmissionAmplitudeByRegion(region,:) = sqrt(layerStack.exitN/layerStack.incidentN)*transmissionAmplitude ...
            .*exp(2i*pi*layerStack.incidentN*referencePlaneGap./wavelengths);
    end
    assert(all(isfinite(transmissionAmplitudeByRegion),'all') && all(abs(transmissionAmplitudeByRegion)<=1+1e-9,'all'), ...
        'UDCImagingV2:FilmPassivity', '膜堆透射系数须有限且满足无源功率约束。');
end

function index = layerIndex(layer,wavelengths,layerImageDir)
    if isfield(layer,'nkFile') && strlength(layer.nkFile)>0
        nkFile = fullfile(layerImageDir,layer.nkFile);
        assert(isfile(nkFile),'UDCImagingV2:NKFile','找不到材料光谱文件：%s',nkFile);
        lines = strtrim(readlines(nkFile));
        lines(ismissing(lines) | strlength(lines)==0) = [];
        if ~isempty(lines) && startsWith(lower(lines(1)),["wavelength","lambda"],'IgnoreCase',true)
            lines(1) = [];
        end
        number = "([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)";
        separator = "(?:\s*,\s*|\s+)";
        rowPattern = "^"+number+separator+number+separator+number+"(?:\s*,)*$";
        opticalConstants = zeros(numel(lines),3);
        for row = 1:numel(lines)
            tokens = regexp(lines(row),rowPattern,'tokens','once');
            assert(~isempty(tokens),'UDCImagingV2:NKTable', ...
                '材料表中存在非三列数值行：%s（数据行 %d）',nkFile,row);
            opticalConstants(row,:) = str2double(tokens);
        end
        tableWavelengthNm = opticalConstants(:,1);
        refractiveIndex = opticalConstants(:,2);
        extinctionCoefficient = opticalConstants(:,3);
        assert(numel(tableWavelengthNm)>=2 && isreal(opticalConstants) && all(isfinite(opticalConstants),'all') && ...
            all(diff(tableWavelengthNm)>0) && all(tableWavelengthNm>0) && ...
            all(refractiveIndex>0) && all(extinctionCoefficient>=0), ...
            'UDCImagingV2:NKTable','材料表须为三列（波长 nm、n、k），波长严格递增，n>0、k≥0：%s',nkFile);
        wavelengthNm = wavelengths*1e9;
        minTableWavelengthNm = tableWavelengthNm(1);
        maxTableWavelengthNm = tableWavelengthNm(end);
        tolerance = 8*eps(max(maxTableWavelengthNm,max(wavelengthNm)));
        assert(all(wavelengthNm>=minTableWavelengthNm-tolerance & wavelengthNm<=maxTableWavelengthNm+tolerance), ...
            'UDCImagingV2:NKRange','材料表未覆盖全部计算波长：%s',nkFile);
        % 消除 m 与 nm 互换的端点舍入误差，实际超出光谱范围的查询已报错。
        wavelengthNm = min(max(wavelengthNm,minTableWavelengthNm),maxTableWavelengthNm);
        index = interp1(tableWavelengthNm,complex(refractiveIndex,extinctionCoefficient),wavelengthNm,'linear');
        return
    end
    validateattributes(layer.n, {'numeric'}, {'vector','real','finite','positive'});
    validateattributes(layer.k, {'numeric'}, {'vector','real','finite','nonnegative'});
    assert((isscalar(layer.n) || numel(layer.n)==numel(wavelengths)) && ...
        (isscalar(layer.k) || numel(layer.k)==numel(wavelengths)), ...
        'UDCImagingV2:IndexSize','n、k 应为标量或与波长一一对应。');
    index = reshape(layer.n,1,[]) + 1i*reshape(layer.k,1,[]);
end
