function [model, geometry] = prepareCoeLayerStack(cfg, stack, inputDir)
%prepareCoeLayerStack 按声明读取图样，统一构造任意层数的局部复反射系数。
% nk 为有限膜；reflectance 为不透明有效反射边界，不由 R 猜测透射系数。
    arguments
        cfg (1,1) struct
        stack (1,1) struct
        inputDir (1,1) string = ""
    end
    names = string(stack.layerNames);
    levels = stack.layerLevels;
    assert(isrow(names) && ~isempty(names) && all(strlength(names)>0) && ...
        numel(unique(names))==numel(names) && all(arrayfun(@(x) isvarname(char(x)), names)), ...
        'COEReflection:LayerNames', 'layerNames 必须为非空、不重复的有效名称列表。');
    assert(~any(ismember(names,["RegionCode","AreaFraction","LocalHeightM"])), ...
        'COEReflection:LayerNames','图层名称不能与区域表统计字段重名。');
    validateattributes(levels, {'numeric'}, {'row','numel',numel(names),'integer','positive','finite'});
    assert(all(diff(levels)>=0), 'COEReflection:LayerOrder', 'layerLevels 必须按从下到上非递减排列。');
    validateattributes(stack.incidentN, {'numeric'}, {'scalar','real','finite','positive'});
    validateattributes(stack.exitN, {'numeric'}, {'scalar','real','finite','positive'});
    mustBeMember(string(stack.backgroundMode), ["medium","absorber"]);
    mustBeMember(string(stack.geometry.source), ["images","mat"]);
    validateattributes(stack.geometry.threshold, {'numeric'}, {'scalar','real','>=',0,'<=',1});
    shape = cfg.periods([2,1])*cfg.samplesPerPeriod;
    masks = false([shape, numel(names)]);
    nk = complex(zeros(numel(names), numel(cfg.wavelengths)));
    thickness = zeros(1,numel(names));
    direct = false(1,numel(names));
    bareCoefficients = complex(zeros(size(nk)));
    data = struct();
    if stack.geometry.source == "mat"
        data = load(fullfile(inputDir, stack.geometry.matFile));
    end
    for j = 1:numel(names)
        assert(isfield(stack.layers,names(j)), 'COEReflection:LayerDefinition', '缺少层定义：%s',names(j));
        layer = stack.layers.(names(j));
        coverage="image";
        if isfield(layer,'coverage'), coverage=string(layer.coverage); end
        mustBeMember(coverage,["image","full"]);
        if coverage=="full"
            value=true(shape); % 可选程序接口；正式IJP项目使用全白图片以统一输入
        elseif stack.geometry.source == "images"
            file = resolveCoeImageFile(fullfile(inputDir,layer.file));
            value = imread(file);
            if isinteger(value)
                value = double(value)/double(intmax(class(value)));
            else
                value = double(value);
            end
            if ndims(value)==3, value=mean(value,3); end
            value = value > stack.geometry.threshold;
        else
            assert(isfield(data,layer.field), 'COEReflection:StackField', 'MAT 缺少图层字段：%s',layer.field);
            value = data.(layer.field);
            assert(all(value(:)==0 | value(:)==1), 'COEReflection:BinaryMask', '图层必须为二值：%s',names(j));
        end
        validateattributes(value, {'numeric','logical'}, {'real','finite','size',shape});
        validateattributes(layer.invert, {'logical'}, {'scalar'});
        masks(:,:,j) = xor(logical(value),layer.invert);
        method = string(layer.method);
        mustBeMember(method,["nk","reflectance"]);
        direct(j) = method=="reflectance";
        if direct(j)
            validateattributes(layer.thickness, {'numeric'}, {'scalar','real','finite','nonnegative'});
            assert(layer.thickness==0, 'COEReflection:ReflectorThickness', ...
                '直接反射率代表不透明的有效反射边界，thickness 必须为 0；有限膜请选 nk。');
            validateattributes(layer.reflectance, {'numeric'}, {'vector','real','finite','>=',0,'<=',1});
            validateattributes(layer.phaseDeg, {'numeric'}, {'vector','real','finite'});
            assert((isscalar(layer.reflectance)||numel(layer.reflectance)==numel(cfg.wavelengths)) && ...
                (isscalar(layer.phaseDeg)||numel(layer.phaseDeg)==numel(cfg.wavelengths)), ...
                'COEReflection:LayerSpectrum', '反射率和相位应为常数或与波长数对应的向量。');
            bareCoefficients(j,:) = sqrt(layer.reflectance(:).').*exp(1i*deg2rad(layer.phaseDeg(:).'));
        else
            [nk(j,:), thickness(j)] = coeMaterialSpectrum(layer,cfg.wavelengths,inputDir);
        end
        geometry.layerMasks.(names(j)) = masks(:,:,j);
    end
    groups = unique(levels);
    referenceHeight = 0;
    for level = groups
        assert(all(sum(masks(:,:,levels==level),3)<=1,'all'), ...
            'COEReflection:SameLevelOverlap', '同一物理层的图片不可交叠，层级 %d。',level);
        referenceHeight = referenceHeight+max(thickness(levels==level));
    end
    % 仅对实际存在的膜堆编码，不构建 2^层数的全组合查找表。
    [patterns,~,regionIndex] = unique(reshape(masks,[],numel(names)),'rows');
    lookup = complex(zeros(size(patterns,1),numel(cfg.wavelengths)));
    localHeight = patterns*thickness(:);
    area = accumarray(regionIndex,1,[size(patterns,1),1])/prod(shape);
    for region = 1:size(patterns,1)
        present = find(patterns(region,:));
        reflector = present(direct(present));
        baseHeight = 0;
        bare = (stack.incidentN-stack.exitN)/(stack.incidentN+stack.exitN);
        if stack.backgroundMode=="absorber", bare=0; end
        if ~isempty(reflector)
            top = reflector(end); % 不透明反射面遮蔽其下方光学膜堆
            bare = bareCoefficients(top,:);
            baseHeight = sum(thickness(present(present<top)));
            present = present(present>top);
        end
        films = fliplr(present(~direct(present))); % 算法输入为从入射侧到下方
        lookup(region,:) = reflectCoeFilmStack(cfg.wavelengths,stack.incidentN, ...
            nk(films,:),thickness(films),bare,referenceHeight-baseHeight);
    end
    model.enabled = true;
    model.regionMap = reshape(uint32(regionIndex-1),shape);
    model.coefficientLookup = lookup;
    model.referenceHeight = referenceHeight;
    model.config = stack;
    model.regionTable = table((0:size(patterns,1)-1).',area,localHeight, ...
        'VariableNames',{'RegionCode','AreaFraction','LocalHeightM'});
    for j = 1:numel(names), model.regionTable.(names(j))=patterns(:,j); end
    material = repelem(names.',numel(cfg.wavelengths));
    methods = strings(numel(names),1);
    for j = 1:numel(names), methods(j)=stack.layers.(names(j)).method; end
    realN=real(nk); imaginaryN=imag(nk); realN(direct,:)=NaN; imaginaryN(direct,:)=NaN;
    columns={material,repelem(levels(:),numel(cfg.wavelengths)),repelem(methods,numel(cfg.wavelengths)), ...
        repmat(cfg.wavelengths(:)*1e9,numel(names),1),reshape(realN.',[],1),reshape(imaginaryN.',[],1), ...
        repelem(thickness(:),numel(cfg.wavelengths)),reshape(abs(bareCoefficients.').^2,[],1), ...
        reshape(rad2deg(angle(bareCoefficients.')),[],1)};
    % 单个图层时 repelem 返回行向量；统一成列，保证多波长材料表形状一致。
    columns=cellfun(@(x)x(:),columns,'UniformOutput',false);
    model.materialTable = table(columns{:},'VariableNames', ...
        {'Material','Level','Method','WavelengthNm','n','k','ThicknessM','DirectReflectance','DirectPhaseDeg'});
    model.materialTable.DirectReflectance(~repelem(direct(:),numel(cfg.wavelengths)))=NaN;
    model.materialTable.DirectPhaseDeg(~repelem(direct(:),numel(cfg.wavelengths)))=NaN;
    geometry.unionMask = any(masks,3);
    geometry.layerNames = names;
    geometry.layerLevels = levels;
    geometry.gridStep = cfg.pitch/cfg.samplesPerPeriod;
end
