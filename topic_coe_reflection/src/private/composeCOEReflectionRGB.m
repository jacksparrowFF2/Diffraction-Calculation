function [withoutHeightRGB, withHeightRGB, mapping] = composeCOEReflectionRGB(withoutHeight, withHeight, options)
%composeCOEReflectionRGB 用共同曝光及对数/伽马映射合成两种高度相位状态的 RGB 预览。
% 三个输入通道必须依次为 R/G/B 代理响应；不进行色度学颜色重建。
% 两种状态、三个通道共用参考值，不分别归一化，保留通道与状态间的相对强弱。
    arguments
        withoutHeight (:, :, 3) double {mustBeReal, mustBeFinite, mustBeNonnegative}
        withHeight (:, :, 3) double {mustBeReal, mustBeFinite, mustBeNonnegative}
        options (1,1) struct
    end
    assert(isequal(size(withoutHeight), size(withHeight)), ...
        'COEReflection:RGBSize', '两种状态的 RGB 信号尺寸必须相同。');
    validateattributes(options.clipQuantile, {'numeric'}, {'scalar', 'real', '>', 0, '<=', 1});
    validateattributes(options.alpha, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
    validateattributes(options.gamma, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
    pooled = sort([withoutHeight(:); withHeight(:)]);
    reference = pooled(max(1, ceil(options.clipQuantile*numel(pooled))));
    if reference == 0
        reference = max(pooled); % 极稀疏或全黑信号也可以显示
    end
    reference = max(reference, realmin);
    mapping = options;
    mapping.referenceSignal = reference;
    withoutHeightRGB = (log1p(options.alpha*min(withoutHeight/reference, 1))/log1p(options.alpha)).^options.gamma;
    withHeightRGB = (log1p(options.alpha*min(withHeight/reference, 1))/log1p(options.alpha)).^options.gamma;
end
