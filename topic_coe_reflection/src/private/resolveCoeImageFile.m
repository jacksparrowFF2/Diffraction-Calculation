function file = resolveCoeImageFile(file)
%resolveCoeImageFile 支持明确的图片文件名，或项目内不带扩展名的标准图层名。
    [folder, stem, extension] = fileparts(file);
    if strlength(extension) == 0
        candidates = dir(fullfile(folder, stem + ".*"));
        supported = [".png", ".bmp", ".jpg", ".jpeg", ".tif", ".tiff"];
        matches = strings(0, 1);
        for j = 1:numel(candidates)
            [~, candidateStem, candidateExtension] = fileparts(candidates(j).name);
            if ~candidates(j).isdir && strcmpi(candidateStem, stem) && any(strcmpi(candidateExtension, supported))
                matches(end+1, 1) = fullfile(folder, candidates(j).name); %#ok<AGROW>
            end
        end
        assert(~isempty(matches), 'COEReflection:StackImage', '未找到图层图片：%s（支持 PNG/BMP/JPG/JPEG/TIF/TIFF）。', file);
        assert(numel(matches) == 1, 'COEReflection:AmbiguousImage', ...
            '同名图层存在多种格式，请在主脚本明确填写扩展名：%s', file);
        file = matches(1);
    end
    assert(isfile(file), 'COEReflection:StackImage', '图样文件不存在：%s', file);
end
