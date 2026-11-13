function [nk, thickness] = coeMaterialSpectrum(material, wavelengths, inputDir)
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
