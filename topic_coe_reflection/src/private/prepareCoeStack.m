function [model, geometry] = prepareCoeStack(cfg, stack, inputDir)
%prepareCoeStack 读取/生成独立层图样，校验覆盖关系，计算各局部膜堆的光谱复系数。
% 不生成 Ny×Nx×波长数的大数组；regionMap + coefficientLookup 保留紧凑表示。
    arguments
        cfg (1,1) struct
        stack (1,1) struct
        inputDir (1,1) string = ""
    end
    validateattributes(cfg.pitch, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
    validateattributes(cfg.periods, {'numeric'}, {'row', 'numel', 2, 'integer', 'positive'});
    validateattributes(cfg.samplesPerPeriod, {'numeric'}, {'scalar', 'integer', '>=', 4});
    validateattributes(stack.incidentN, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
    validateattributes(stack.anodeReflectance, {'numeric'}, {'scalar', 'real', '>=', 0, '<=', 1});
    validateattributes(stack.anodePhaseDeg, {'numeric'}, {'scalar', 'real', 'finite'});
    validateattributes(stack.bmExpansion, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'});
    shape = [cfg.periods(2), cfg.periods(1)]*cfg.samplesPerPeriod;
    dx = cfg.pitch/cfg.samplesPerPeriod;
    geometry = readGeometry(cfg, stack, inputDir, shape, dx);
    assert(any(geometry.anodeMask, 'all'), 'COEReflection:EmptyAnode', '阳极图样不能为空。');
    assert(all(geometry.bmOpenMask(geometry.anodeMask)), ...
        'COEReflection:BMContainment', 'BM 开口必须包含阳极图样。');
    cfCount = double(geometry.cfRMask) + double(geometry.cfGMask) + double(geometry.cfBMask);
    assert(all(cfCount <= 1, 'all'), 'COEReflection:CFOverlap', 'RGB CF 图样存在交叠，请缩小填充或校正对齐。');
    if stack.geometry.requireCFCoverage
        assert(all(cfCount(geometry.bmOpenMask) == 1), ...
            'COEReflection:CFCoverage', 'RGB CF 必须完整覆盖 BM 开口。');
        for name = ["cfRMask", "cfGMask", "cfBMask"]
            mask = geometry.(name);
            assert(~any(mask, 'all') || any(mask & ~geometry.bmOpenMask, 'all'), ...
                'COEReflection:CFExtent', '非空 %s 应向 BM 实体上方延伸，不能仅等于开口。', name);
        end
    end
    geometry.cfLabel = uint8(geometry.cfRMask) + 2*uint8(geometry.cfGMask) + 3*uint8(geometry.cfBMask);
    geometry.absorberMask = ~geometry.anodeMask;
    geometry.expansionAbsorberMask = geometry.bmOpenMask & ~geometry.anodeMask;
    geometry.bmMask = ~geometry.bmOpenMask;
    regionMap = uint8(geometry.anodeMask) + 2*uint8(geometry.bmMask) + 4*geometry.cfLabel;
    regionIDs = unique(regionMap(:));
    names = ["BM", "R", "G", "B"];
    nk = zeros(4, numel(cfg.wavelengths));
    depths = zeros(1, 4);
    for j = 1:4
        [nk(j, :), depths(j)] = materialSpectrum(stack.materials.(names(j)), cfg.wavelengths, inputDir);
    end
    referenceHeight = depths(1) + max(depths(2:4));
    lookup = complex(zeros(16, numel(cfg.wavelengths)));
    area = zeros(numel(regionIDs), 1);
    hasAnode = false(size(area)); hasBM = hasAnode; cf = zeros(size(area));
    for j = 1:numel(regionIDs)
        code = double(regionIDs(j));
        hasAnode(j) = bitand(code, 1) ~= 0;
        hasBM(j) = bitand(code, 2) ~= 0;
        cf(j) = floor(code/4);
        layerIDs = [cf(j) + 1, 1];
        layerIDs = layerIDs([cf(j) > 0, hasBM(j)]);
        bare = 0; % 阳极外：裸底零反射的匹配吸收终端，而不是扩展阳极
        if hasAnode(j)
            bare = sqrt(stack.anodeReflectance)*exp(1i*deg2rad(stack.anodePhaseDeg));
        end
        lookup(code + 1, :) = reflectCoeFilmStack(cfg.wavelengths, stack.incidentN, ...
            nk(layerIDs, :), depths(layerIDs), bare, referenceHeight);
        area(j) = nnz(regionMap == code)/numel(regionMap);
    end
    model.enabled = true;
    model.regionMap = regionMap;
    model.coefficientLookup = lookup;
    model.referenceHeight = referenceHeight;
    model.config = stack;
    model.regionTable = table(double(regionIDs), hasAnode, hasBM, cf, area, ...
        'VariableNames', {'RegionCode', 'HasAnode', 'HasBM', 'CFLabel', 'AreaFraction'});
    model.materialTable = table(repelem(names.', numel(cfg.wavelengths)), ...
        repmat(cfg.wavelengths(:)*1e9, 4, 1), reshape(real(nk).', [], 1), ...
        reshape(imag(nk).', [], 1), repelem(depths.', numel(cfg.wavelengths)), ...
        'VariableNames', {'Material', 'WavelengthNm', 'n', 'k', 'ThicknessM'});
    geometry.gridStep = dx;
end

function geometry = readGeometry(cfg, stack, inputDir, shape, dx)
    source = string(stack.geometry.source);
    mustBeMember(source, ["diamond", "mat", "images"]);
    switch source
        case "diamond"
            geometry = createCoeDiamondGeometry(cfg, stack);
        case "mat"
            file = fullfile(inputDir, stack.geometry.matFile);
            assert(isfile(file), 'COEReflection:StackInput', '多层图样 MAT 文件不存在：%s', file);
            available = string(who('-file', file));
            requested = ["anodeMask", "bmOpenMask", "bmMask", "cfRMask", "cfGMask", "cfBMask"];
            selected = cellstr(intersect(requested, available, 'stable'));
            data = struct();
            if ~isempty(selected)
                data = load(file, selected{:});
            end
            fields = ["anodeMask", "cfRMask", "cfGMask", "cfBMask"];
            for name = fields
                assert(isfield(data, name), 'COEReflection:StackField', 'MAT 文件缺少字段 %s。', name);
                geometry.(name) = binaryMask(data.(name), shape, name);
            end
            if isfield(data, 'bmOpenMask')
                geometry.bmOpenMask = binaryMask(data.bmOpenMask, shape, 'bmOpenMask');
                if isfield(data, 'bmMask')
                    assert(isequal(~geometry.bmOpenMask, binaryMask(data.bmMask, shape, 'bmMask')), ...
                        'COEReflection:BMAmbiguous', 'bmOpenMask 与 bmMask 必须互为补集。');
                end
            elseif isfield(data, 'bmMask')
                geometry.bmOpenMask = ~binaryMask(data.bmMask, shape, 'bmMask');
            else
                geometry.bmOpenMask = expandMask(geometry.anodeMask, stack.bmExpansion, dx);
            end
        case "images"
            geometry.anodeMask = imageMask(fullfile(inputDir, stack.geometry.anodeFile), shape, stack.geometry.threshold);
            geometry.cfRMask = imageMask(fullfile(inputDir, stack.geometry.cfFiles(1)), shape, stack.geometry.threshold);
            geometry.cfGMask = imageMask(fullfile(inputDir, stack.geometry.cfFiles(2)), shape, stack.geometry.threshold);
            geometry.cfBMask = imageMask(fullfile(inputDir, stack.geometry.cfFiles(3)), shape, stack.geometry.threshold);
            if strlength(stack.geometry.bmFile) == 0
                geometry.bmOpenMask = expandMask(geometry.anodeMask, stack.bmExpansion, dx);
            else
                geometry.bmOpenMask = imageMask(fullfile(inputDir, stack.geometry.bmFile), shape, stack.geometry.threshold);
                if ~stack.geometry.bmImageIsOpening
                    geometry.bmOpenMask = ~geometry.bmOpenMask;
                end
            end
    end
end

function mask = binaryMask(value, shape, name)
    validateattributes(value, {'numeric', 'logical'}, {'real', 'finite', 'size', shape}, '', name);
    assert(all(value(:) == 0 | value(:) == 1), 'COEReflection:BinaryMask', '%s 必须为二值 0/1。', name);
    mask = logical(value);
end

function mask = imageMask(file, shape, threshold)
    file = resolveCoeImageFile(file);
    validateattributes(threshold, {'numeric'}, {'scalar', 'real', '>=', 0, '<=', 1});
    value = imread(file);
    if isinteger(value)
        value = double(value)/double(intmax(class(value)));
    else
        value = double(value);
    end
    if ndims(value) == 3
        value = mean(value, 3);
    end
    mask = binaryMask(value > threshold, shape, 'imageMask');
end

function expanded = expandMask(mask, distance, dx)
    radius = ceil(distance/dx);
    [xx, yy] = meshgrid(-radius:radius);
    disk = hypot(xx, yy)*dx <= distance + 1e-15;
    expanded = conv2(double(mask), double(disk), 'same') > 0;
end

function [nk, thickness] = materialSpectrum(material, wavelengths, inputDir)
    validateattributes(material.thickness, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'});
    thickness = material.thickness;
    if isfield(material, 'nkFile') && strlength(material.nkFile) > 0
        data = readmatrix(fullfile(inputDir, material.nkFile));
        assert(size(data, 2) == 3 && all(isfinite(data), 'all') && all(diff(data(:, 1)) > 0), ...
            'COEReflection:NKTable', 'n/k 文件必须为严格递增的三列：波长 nm、n、k。');
        n = interp1(data(:, 1), data(:, 2), wavelengths*1e9, 'linear');
        k = interp1(data(:, 1), data(:, 3), wavelengths*1e9, 'linear');
    else
        n = material.n; k = material.k;
    end
    validateattributes(n, {'numeric'}, {'vector', 'real', 'finite', 'positive'});
    validateattributes(k, {'numeric'}, {'vector', 'real', 'finite', 'nonnegative'});
    assert((isscalar(n) || numel(n) == numel(wavelengths)) && (isscalar(k) || numel(k) == numel(wavelengths)), ...
        'COEReflection:NKSamples', 'n、k 应为常数或与波长一一对应的向量，文件须覆盖全部计算波段。');
    nk = reshape(n, 1, []).*ones(size(wavelengths)) + 1i*reshape(k, 1, []).*ones(size(wavelengths));
end
