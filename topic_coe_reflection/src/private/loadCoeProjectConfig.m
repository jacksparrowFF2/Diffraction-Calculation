function [settings, record] = loadCoeProjectConfig(file, settings)
%loadCoeProjectConfig 读取 UTF-8 key=value 文本，仅允许覆盖已有配置字段。
% 数字/数组/布尔/字符串使用 JSON 字面量；不执行 MATLAB 表达式或脚本。
    assert(isfile(file), 'COEReflection:MissingConfig', '项目配置不存在：%s', file);
    record.file = string(file);
    record.text = string(fileread(file));
    lines = splitlines(record.text);
    % 先收集声明列表，让自定义层名的字段不受文本先后顺序影响；仍严格校验字段名。
    if isfield(settings,'stack') && isfield(settings.stack,'layerNames')
        for index = 1:numel(lines)
            declaration = regexp(char(strtrim(lines(index))), '^stack\.layerNames\s*=\s*(.+)$','tokens','once');
            if ~isempty(declaration)
                names = reshape(string(jsondecode(declaration{1})),1,[]);
                assert(~isempty(names) && numel(unique(names))==numel(names) && ...
                    all(arrayfun(@(x) isvarname(char(x)),names)), ...
                    'COEReflection:LayerNames','layerNames 应为非空、唯一的有效名称列表。');
                for name = names
                    if ~isfield(settings.stack.layers,name)
                        settings.stack.layers.(name)=coeLayerDefaults(name);
                    end
                end
                settings.stack.layerNames=names;
                settings.stack.layerLevels=1:numel(names);
            end
        end
    end
    record.keys = strings(0, 1);
    for index = 1:numel(lines)
        line = strtrim(lines(index));
        if index == 1
            line = erase(line, char(65279)); % 接受带 BOM 的 UTF-8 文本
        end
        if strlength(line) == 0 || startsWith(line, ["#", "%"])
            continue
        end
        tokens = regexp(char(line), '^([A-Za-z]\w*(?:\.[A-Za-z]\w*)*)\s*=\s*(.+)$', 'tokens', 'once');
        assert(~isempty(tokens), 'COEReflection:ConfigSyntax', '配置第 %d 行应为 key=value：%s', index, line);
        key = string(tokens{1});
        assert(~any(key == ["cfg.radius", "cfg.aspectRatio", "cfg.cellAnglesDeg"]), ...
            'COEReflection:UnknownConfig', '像素轮廓请通过图片输入，不支持形状配置：%s', key);
        assert(~any(record.keys == key), 'COEReflection:DuplicateConfig', '配置字段重复：%s', key);
        try
            value = jsondecode(tokens{2});
        catch
            error('COEReflection:ConfigValue', '配置第 %d 行的值无效，请使用数字、JSON 数组、true/false 或双引号字符串。', index);
        end
        if iscell(value) && all(cellfun(@ischar, value(:)))
            value = reshape(string(value), 1, []);
        end
        assert(~isstruct(value) && ~iscell(value), 'COEReflection:ConfigValue', '不支持对象或混合数组：%s', key);
        if ischar(value) || isstring(value)
            value = string(value);
        elseif isnumeric(value) && isvector(value)
            value = reshape(value, 1, []); % JSON 向量统一为行向量，与主脚本约定一致
        end
        if key == "stack.geometry.source"
            assert(isstring(value) && isscalar(value) && any(value == ["images", "mat"]), ...
                'COEReflection:ConfigValue', '正式项目仅支持 images 或 mat 图样输入。');
        end
        settings = replaceField(settings, split(key, '.'), value, key);
        record.keys(end+1, 1) = key;
    end
end

function settings = replaceField(settings, names, value, key)
    name = char(names(1));
    assert(isstruct(settings) && isscalar(settings) && isfield(settings, name), ...
        'COEReflection:UnknownConfig', '未知配置字段：%s', key);
    if numel(names) > 1
        settings.(name) = replaceField(settings.(name), names(2:end), value, key);
    else
        original = settings.(name);
        valid = (isnumeric(original) && isnumeric(value)) || ...
            (islogical(original) && islogical(value)) || ...
            ((isstring(original) || ischar(original)) && isstring(value));
        assert(valid && ~isstruct(original), 'COEReflection:ConfigType', '配置字段类型不匹配：%s', key);
        if ischar(original)
            value = char(value);
        end
        settings.(name) = value;
    end
end
