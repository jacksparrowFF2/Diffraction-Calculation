function [model, geometry] = corePrepareLayerStack(cfg, stack, inputDir)
%corePrepareLayerStack 按声明读取图样，统一构造任意层数的局部复反射系数。
% nk 为有限膜；reflectance 为不透明有效反射边界，不由 R 猜测透射系数。
    arguments
        cfg (1,1) struct
        stack (1,1) struct
        inputDir (1,1) string = ""
    end
    [geometry, regionMap, patterns] = displayoptics.readLayerGeometry(cfg, stack, inputDir);
    names = geometry.layerNames;
    levels = geometry.layerLevels;
    validateattributes(stack.incidentN, {'numeric'}, {'scalar','real','finite','positive'});
    validateattributes(stack.exitN, {'numeric'}, {'scalar','real','finite','positive'});
    mustBeMember(string(stack.backgroundMode), ["medium","absorber"]);
    nk = complex(zeros(numel(names), numel(cfg.wavelengths)));
    thickness = zeros(1,numel(names));
    direct = false(1,numel(names));
    bareCoefficients = complex(zeros(size(nk)));
    for j = 1:numel(names)
        layer = stack.layers.(names(j));
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
            % 材料光谱优先读取 nkFile，否则使用常数或逐波长 n/k。
            validateattributes(layer.thickness, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'});
            thickness(j) = layer.thickness;
            if isfield(layer, 'nkFile') && strlength(layer.nkFile) > 0
                nkData = readmatrix(fullfile(inputDir, layer.nkFile));
                assert(size(nkData, 2) == 3 && all(isfinite(nkData), 'all') && all(diff(nkData(:, 1)) > 0), ...
                    'COEReflection:NKTable', 'n/k 文件必须为严格递增的三列：波长 nm、n、k。');
                n = interp1(nkData(:, 1), nkData(:, 2), cfg.wavelengths*1e9, 'linear');
                k = interp1(nkData(:, 1), nkData(:, 3), cfg.wavelengths*1e9, 'linear');
            else
                n = layer.n; k = layer.k;
            end
            validateattributes(n, {'numeric'}, {'vector', 'real', 'finite', 'positive'});
            validateattributes(k, {'numeric'}, {'vector', 'real', 'finite', 'nonnegative'});
            assert((isscalar(n) || numel(n) == numel(cfg.wavelengths)) && ...
                (isscalar(k) || numel(k) == numel(cfg.wavelengths)), ...
                'COEReflection:NKSamples', 'n、k 应为常数或与波长一一对应的向量，文件须覆盖全部计算波段。');
            nk(j,:) = reshape(n, 1, []).*ones(size(cfg.wavelengths)) ...
                + 1i*reshape(k, 1, []).*ones(size(cfg.wavelengths));
        end
    end
    referenceHeight = 0;
    for level = unique(levels)
        referenceHeight = referenceHeight+max(thickness(levels==level));
    end
    lookup = complex(zeros(size(patterns,1),numel(cfg.wavelengths)));
    localHeight = patterns*thickness(:);
    area = accumarray(double(regionMap(:))+1,1,[size(patterns,1),1])/numel(regionMap);
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
        lookup(region,:) = coreReflectFilmStack(cfg.wavelengths,stack.incidentN, ...
            nk(films,:),thickness(films),bare,referenceHeight-baseHeight);
    end
    model.enabled = true;
    model.regionMap = regionMap;
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
end
