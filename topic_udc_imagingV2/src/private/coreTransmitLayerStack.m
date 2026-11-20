function lookup = coreTransmitLayerStack(wavelengths, stack, patterns, inputDir)
%coreTransmitLayerStack 求各区域法向复透射，以入射功率归一化振幅。
    validateattributes(stack.incidentN, {'numeric'}, {'scalar','real','finite','positive'});
    validateattributes(stack.exitN, {'numeric'}, {'scalar','real','finite','positive'});
    names = string(stack.layerNames);
    background = "opaque";
    if isfield(stack,'backgroundMode'), background = string(stack.backgroundMode); end
    mustBeMember(background,["opaque","medium"]);
    direct = false(1,numel(names));
    directAmplitude = zeros(numel(names),numel(wavelengths));
    nk = zeros(numel(names),numel(wavelengths));
    thickness = zeros(1,numel(names));
    for j = 1:numel(names)
        layer = stack.layers.(names(j));
        mustBeMember(string(layer.method), ["transmittance","nk"]);
        direct(j) = string(layer.method) == "transmittance";
        if direct(j)
            validateattributes(layer.transmittance, {'numeric'}, {'vector','real','finite','>=',0,'<=',1});
            assert(isscalar(layer.transmittance) || numel(layer.transmittance)==numel(wavelengths), ...
                'UDCImagingV2:TransmittanceSize', '透过率应为标量或与波长一一对应。');
            directAmplitude(j,:) = sqrt(reshape(layer.transmittance,1,[]));
        else
            validateattributes(layer.thickness, {'numeric'}, {'scalar','real','finite','nonnegative'});
            nk(j,:) = layerIndex(layer,wavelengths,inputDir);
            thickness(j) = layer.thickness;
        end
    end
    referenceHeight = 0;
    for level = unique(stack.layerLevels)
        referenceHeight = referenceHeight + max(thickness(stack.layerLevels==level));
    end

    lookup = complex(zeros(size(patterns,1),numel(wavelengths)));
    for region = 1:size(patterns,1)
        present = find(patterns(region,:));
        if isempty(present) && background=="opaque"
            continue
        end
        effective = present(direct(present));
        if ~isempty(effective)
            assert(isscalar(present), 'UDCImagingV2:DirectOverlap', ...
                '直接透过率代表区域的整体透射，白区不可与其他透射层重叠。');
            lookup(region,:) = directAmplitude(effective,:);
            continue
        end
        % 图层从下到上存储，从出射界面向入射侧递推。
        films = find(patterns(region,:) & thickness>0);
        indices = [stack.exitN*ones(size(wavelengths)); nk(films,:); stack.incidentN*ones(size(wavelengths))];
        r = zeros(size(wavelengths));
        t = ones(size(wavelengths));
        p = t;
        for j = 2:size(indices,1)
            lower = indices(j-1,:);
            upper = indices(j,:);
            interface = (upper-lower)./(upper+lower);
            denominator = 1 + interface.*r.*p.^2;
            t = (1+interface).*p.*t./denominator;
            r = (interface+r.*p.^2)./denominator;
            if j <= numel(films)+1
                p = exp(2i*pi*upper*thickness(films(j-1))./wavelengths);
            end
        end
        gap = max(0,referenceHeight-sum(thickness(films)));
        lookup(region,:) = sqrt(stack.exitN/stack.incidentN)*t ...
            .*exp(2i*pi*stack.incidentN*gap./wavelengths);
    end
    assert(all(isfinite(lookup),'all') && all(abs(lookup)<=1+1e-9,'all'), ...
        'UDCImagingV2:FilmPassivity', '膜堆透射系数须有限且满足无源功率约束。');
end

function index = layerIndex(layer,wavelengths,inputDir)
    if isfield(layer,'nkFile') && strlength(layer.nkFile)>0
        file = fullfile(inputDir,layer.nkFile);
        assert(isfile(file),'UDCImagingV2:NKFile','找不到材料光谱文件：%s',file);
        data = readNkRows(file);
        assert(size(data,2)==3 && size(data,1)>=2 && isreal(data) && all(isfinite(data),'all') && ...
            all(diff(data(:,1))>0) && all(data(:,1)>0) && all(data(:,2)>0) && all(data(:,3)>=0), ...
            'UDCImagingV2:NKTable','材料表须为三列（波长 nm、n、k），波长严格递增，n>0、k≥0：%s',file);
        wavelengthNm = wavelengths*1e9;
        lower = data(1,1);
        upper = data(end,1);
        tolerance = 8*eps(max(upper,max(wavelengthNm)));
        assert(all(wavelengthNm>=lower-tolerance & wavelengthNm<=upper+tolerance), ...
            'UDCImagingV2:NKRange','材料表未覆盖全部计算波长：%s',file);
        % 消除 m 与 nm 互换的端点舍入误差，实际超出光谱范围的查询已报错。
        wavelengthNm = min(max(wavelengthNm,lower),upper);
        index = interp1(data(:,1),complex(data(:,2),data(:,3)),wavelengthNm,'linear');
        return
    end
    validateattributes(layer.n, {'numeric'}, {'vector','real','finite','positive'});
    validateattributes(layer.k, {'numeric'}, {'vector','real','finite','nonnegative'});
    assert((isscalar(layer.n) || numel(layer.n)==numel(wavelengths)) && ...
        (isscalar(layer.k) || numel(layer.k)==numel(wavelengths)), ...
        'UDCImagingV2:IndexSize','n、k 应为标量或与波长一一对应。');
    index = reshape(layer.n,1,[]) + 1i*reshape(layer.k,1,[]);
end

function data = readNkRows(file)
    lines = strtrim(readlines(file));
    lines(ismissing(lines) | strlength(lines)==0) = [];
    if isempty(lines), data = zeros(0,3); return; end
    if startsWith(lower(lines(1)),["wavelength","lambda"],'IgnoreCase',true)
        lines(1) = [];
    end
    number = "([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)";
    separator = "(?:\s*,\s*|\s+)";
    pattern = "^"+number+separator+number+separator+number+"(?:\s*,)*$";
    data = zeros(numel(lines),3);
    for row = 1:numel(lines)
        tokens = regexp(lines(row),pattern,'tokens','once');
        assert(~isempty(tokens),'UDCImagingV2:NKTable','材料表中存在非三列数值行：%s（数据行 %d）',file,row);
        data(row,:) = str2double(tokens);
    end
end
