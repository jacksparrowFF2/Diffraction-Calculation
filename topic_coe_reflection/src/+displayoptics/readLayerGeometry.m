function [geometry, regionMap, patterns] = readLayerGeometry(cfg, stack, inputDir, strictBinary)
%readLayerGeometry 读取并对齐图层，按实际覆盖组合编码物面区域。
    if nargin < 4, strictBinary = false; end
    names = string(stack.layerNames);
    levels = stack.layerLevels;
    assert(isrow(names) && ~isempty(names) && all(strlength(names)>0) && ...
        numel(unique(names))==numel(names) && all(arrayfun(@(x) isvarname(char(x)), names)), ...
        'COEReflection:LayerNames', 'layerNames 必须为非空、不重复的有效名称列表。');
    assert(~any(ismember(names,["RegionCode","AreaFraction","LocalHeightM"])), ...
        'COEReflection:LayerNames','图层名称不能与区域表统计字段重名。');
    validateattributes(levels, {'numeric'}, {'row','numel',numel(names),'integer','positive','finite'});
    assert(all(diff(levels)>=0), 'COEReflection:LayerOrder', 'layerLevels 必须按从下到上非递减排列。');
    mustBeMember(string(stack.geometry.source), ["images","mat"]);
    if ~strictBinary
        validateattributes(stack.geometry.threshold, {'numeric'}, {'scalar','real','>=',0,'<=',1});
    end
    shape = cfg.periods([2,1])*cfg.samplesPerPeriod;
    repeatToFit=false;
    if isfield(stack.geometry,'repeatToFit'), repeatToFit=stack.geometry.repeatToFit; end
    validateattributes(repeatToFit,{'logical'},{'scalar'});
    cropToFit=false; pixelReplication=false;
    if isfield(stack.geometry,'cropToFit'), cropToFit=stack.geometry.cropToFit; end
    if isfield(stack.geometry,'pixelReplication'), pixelReplication=stack.geometry.pixelReplication; end
    validateattributes(cropToFit,{'logical'},{'scalar'});
    validateattributes(pixelReplication,{'logical'},{'scalar'});
    assert((~cropToFit && ~pixelReplication)||stack.geometry.source=="images", ...
        'COEReflection:RasterSource','显式裁取/像素复制仅适用于图片输入。');
    inputSamples=cfg.samplesPerPeriod;
    if isfield(stack.geometry,'inputSamplesPerPeriod')
        validateattributes(stack.geometry.inputSamplesPerPeriod,{'numeric'}, ...
            {'scalar','integer','nonnegative','finite'});
        if stack.geometry.inputSamplesPerPeriod>0, inputSamples=stack.geometry.inputSamplesPerPeriod; end
    end
    validateattributes(inputSamples,{'numeric'},{'scalar','integer','positive','finite'});
    assert(stack.geometry.source=="images" || inputSamples==cfg.samplesPerPeriod, ...
        'COEReflection:RasterSource','输入采样率重采样仅适用于图片输入。');
    geometry.inputSamplesPerPeriod=inputSamples;
    geometry.samplesPerPeriod=cfg.samplesPerPeriod;
    assert(~repeatToFit || stack.geometry.source=="images", ...
        'COEReflection:RepeatSource','repeatToFit只用于图片输入。');
    inputShape=[];
    geometry.inputPeriods=cfg.periods;
    geometry.repeatCounts=[1,1];
    geometry.inputGridSize=shape;
    geometry.cropBoundsPixels=[1,shape(1),1,shape(2)];
    geometry.selectedPeriods=cfg.periods;
    geometry.pixelReplicationFactor=1;
    masks = false([shape, numel(names)]);
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
            % 支持明确文件名或无扩展名图层名；同名多格式必须显式选择。
            file = fullfile(inputDir,layer.file);
            [folder, stem, extension] = fileparts(file);
            if strlength(extension) == 0
                candidates = dir(fullfile(folder, stem + ".*"));
                supported = [".png", ".bmp", ".jpg", ".jpeg", ".tif", ".tiff"];
                matches = strings(0, 1);
                for candidateIndex = 1:numel(candidates)
                    candidate = candidates(candidateIndex);
                    [~, candidateStem, candidateExtension] = fileparts(candidate.name);
                    if ~candidate.isdir && strcmpi(candidateStem, stem) && any(strcmpi(candidateExtension, supported))
                        matches(end+1, 1) = fullfile(folder, candidate.name); %#ok<AGROW>
                    end
                end
                assert(~isempty(matches), 'COEReflection:StackImage', ...
                    '未找到图层图片：%s（支持 PNG/BMP/JPG/JPEG/TIF/TIFF）。', file);
                assert(isscalar(matches), 'COEReflection:AmbiguousImage', ...
                    '同名图层存在多种格式，请在主脚本明确填写扩展名：%s', file);
                file = matches(1);
            end
            assert(isfile(file), 'COEReflection:StackImage', '图样文件不存在：%s', file);
            if strictBinary
                value = readBinaryImage(file);
            else
                value = imread(file);
                if isinteger(value)
                    value = double(value)/double(intmax(class(value)));
                else
                    value = double(value);
                end
                if ndims(value)==3, value=mean(value,3); end
                value = value > stack.geometry.threshold;
            end
        else
            assert(isfield(data,layer.field), 'COEReflection:StackField', 'MAT 缺少图层字段：%s',layer.field);
            value = data.(layer.field);
            assert(all(value(:)==0 | value(:)==1), 'COEReflection:BinaryMask', '图层必须为二值：%s',names(j));
        end
        if coverage=="image"
            actualShape=size(value);
            if isempty(inputShape)
                inputShape=actualShape;
                geometry.inputGridSize=inputShape;
                geometry.inputPeriods=inputShape([2,1])/inputSamples;
            end
            assert(isequal(actualShape,inputShape),'COEReflection:InputAlignment', ...
                '所有图层原始图片必须同尺寸、同原点，不能各自采用不同单元。');
            if inputSamples~=cfg.samplesPerPeriod
                assert(all(mod(inputShape,inputSamples)==0), ...
                    'COEReflection:InputPitch','原始图片行列数须为inputSamplesPerPeriod的整数倍。');
                assert(cfg.samplesPerPeriod<inputSamples || ...
                    (pixelReplication && mod(cfg.samplesPerPeriod,inputSamples)==0),'COEReflection:RasterUpsampling', ...
                    '默认禁止提高图片采样率；显式pixelReplication仅允许整数倍复制，不增加真实版图细节。');
                sampledShape=inputShape/inputSamples*cfg.samplesPerPeriod;
                rows=floor(((0:sampledShape(1)-1)+0.5)*inputShape(1)/sampledShape(1))+1;
                cols=floor(((0:sampledShape(2)-1)+0.5)*inputShape(2)/sampledShape(2))+1;
                value=value(rows,cols); % 共同像素中心最近邻采样，保持二值区域及坐标尺度。
                if cfg.samplesPerPeriod>inputSamples
                    geometry.pixelReplicationFactor=cfg.samplesPerPeriod/inputSamples;
                end
            end
            sampledShape=size(value);
            geometry.sampledInputGridSize=sampledShape;
            selectedShape=sampledShape;
            if cropToFit
                selectedShape=min(shape,sampledShape);
                first=(sampledShape-selectedShape)/2+1;
                assert(all(mod(first-1,cfg.samplesPerPeriod)==0) && all(mod(selectedShape,cfg.samplesPerPeriod)==0), ...
                    'COEReflection:CropPitch','中心裁取必须沿完整Pitch边界，不能改变像素尺寸。');
                last=first+selectedShape-1;
                value=value(first(1):last(1),first(2):last(2));
                geometry.cropBoundsPixels=[first(1),last(1),first(2),last(2)];
            else
                geometry.cropBoundsPixels=[1,sampledShape(1),1,sampledShape(2)];
            end
            geometry.selectedPeriods=selectedShape([2,1])/cfg.samplesPerPeriod;
            if repeatToFit
                assert(all(mod(sampledShape,cfg.samplesPerPeriod)==0), ...
                    'COEReflection:InputPitch','图片行列数须为samplesPerPeriod的整数倍，才能确定输入Pitch数。');
                assert(all(mod(shape,selectedShape)==0), ...
                    'COEReflection:RepeatExtent','目标Pitch数须为输入单元Pitch数的整数倍，不缩放或裁切输入。');
                repeats=shape./selectedShape;
                geometry.repeatCounts=repeats([2,1]);
                value=repmat(value,repeats);
            end
        end
        validateattributes(value, {'numeric','logical'}, {'real','finite','size',shape});
        validateattributes(layer.invert, {'logical'}, {'scalar'});
        assert(~strictBinary || ~layer.invert, ...
            'UDCImagingV2:MaskPolarity', '透射图层的 invert 须为 false，请用 materialOn 指定材料覆盖颜色：%s', names(j));
        selectBlack = layer.invert;
        if strictBinary && isfield(layer,'materialOn')
            color = string(layer.materialOn);
            mustBeMember(color,["black","white"]);
            selectBlack = color=="black";
        end
        masks(:,:,j) = xor(logical(value),selectBlack);
        geometry.layerMasks.(names(j)) = masks(:,:,j);
    end
    for level = unique(levels)
        assert(all(sum(masks(:,:,levels==level),3)<=1,'all'), ...
            'COEReflection:SameLevelOverlap', '同一物理层的图片不可交叠，层级 %d。',level);
    end
    % 仅对实际存在的膜堆编码，不构建 2^层数的全组合查找表。
    [patterns,~,regionIndex] = unique(reshape(masks,[],numel(names)),'rows');
    regionMap = reshape(uint32(regionIndex-1),shape);
    geometry.unionMask = any(masks,3);
    geometry.layerNames = names;
    geometry.layerLevels = levels;
    geometry.gridStep = cfg.pitch/cfg.samplesPerPeriod;
end

function value = readBinaryImage(file)
    [image, palette] = imread(file);
    if ~isempty(palette)
        image = ind2rgb(image,palette);
    end
    white = 1;
    if isinteger(image), white = intmax(class(image)); end
    blackPixels = all(image==0,3);
    whitePixels = all(image==white,3);
    assert(all(blackPixels | whitePixels,'all'), ...
        'UDCImagingV2:BinaryMask', '透射掩膜只接受纯黑（0）和纯白（1）：%s', file);
    value = whitePixels;
end
