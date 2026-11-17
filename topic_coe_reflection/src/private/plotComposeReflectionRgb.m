function [withoutHeightRGB, withHeightRGB, mapping] = plotComposeReflectionRgb(withoutHeight, withHeight, options, compareStates)
%plotComposeReflectionRgb 用共同曝光及对数/伽马映射合成两种高度相位状态的 RGB 预览。
% 三个输入通道必须依次为 R/G/B 代理响应；不进行色度学颜色重建。
% 两种状态、三个通道共用参考值，不分别归一化，保留通道与状态间的相对强弱。
    arguments
        withoutHeight (:, :, 3) double {mustBeReal, mustBeFinite, mustBeNonnegative}
        withHeight (:, :, 3) double {mustBeReal, mustBeFinite, mustBeNonnegative}
        options (1,1) struct
        compareStates (1,1) logical = true
    end
    assert(isequal(size(withoutHeight), size(withHeight)), ...
        'COEReflection:RGBSize', '两种状态的 RGB 信号尺寸必须相同。');
    validateattributes(options.clipQuantile, {'numeric'}, {'scalar', 'real', '>', 0, '<=', 1});
    validateattributes(options.alpha, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
    validateattributes(options.gamma, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
    if isfield(options,'fixedReferenceSignal')
        validateattributes(options.fixedReferenceSignal,{'numeric'},{'scalar','real','finite','nonnegative'});
    end
    if ~isfield(options,'mode'), options.mode="legacy"; end
    if ~isfield(options,'peakLevel'), options.peakLevel=1; end
    mustBeMember(options.mode,["legacy","luminance"]);
    validateattributes(options.peakLevel,{'numeric'},{'scalar','real','finite','>',0,'<=',1});
    if options.mode=="luminance"
        withoutY=sum(withoutHeight.*reshape([0.2126,0.7152,0.0722],1,1,3),3);
        withY=sum(withHeight.*reshape([0.2126,0.7152,0.0722],1,1,3),3);
        pooled=withoutY(:);
        if compareStates, pooled=[pooled;withY(:)]; end
        pooled=sort(pooled);
        reference=pooled(max(1,ceil(options.clipQuantile*numel(pooled))));
        if reference==0, reference=max(pooled); end
        reference=max(reference,realmin);
        if isfield(options,'fixedReferenceSignal') && options.fixedReferenceSignal>0
            reference=options.fixedReferenceSignal;
        end
        mapping=options; mapping.referenceSignal=reference;
        mapping.method="共同亮度参考；单调对数映射和RGB向量统一高光压缩；无位置遮罩，原始能量不变";
        withoutHeightRGB=renderLuminance(withoutHeight,withoutY,reference,options);
        if compareStates
            withHeightRGB=renderLuminance(withHeight,withY,reference,options);
        else
            withHeightRGB=withoutHeightRGB;
        end
        return
    end
    if compareStates
        pooled = sort([withoutHeight(:); withHeight(:)]);
    else
        pooled = sort(withoutHeight(:));
    end
    reference = pooled(max(1, ceil(options.clipQuantile*numel(pooled))));
    if reference == 0
        reference = max(pooled); % 极稀疏或全黑信号也可以显示
    end
    reference = max(reference, realmin);
    if isfield(options,'fixedReferenceSignal') && options.fixedReferenceSignal>0
        reference=options.fixedReferenceSignal;
    end
    mapping = options;
    mapping.referenceSignal = reference;
    withoutHeightRGB = (log1p(options.alpha*min(withoutHeight/reference, 1))/log1p(options.alpha)).^options.gamma;
    if compareStates
        withHeightRGB = (log1p(options.alpha*min(withHeight/reference, 1))/log1p(options.alpha)).^options.gamma;
    else
        withHeightRGB = withoutHeightRGB;
    end
end

function rgb=renderLuminance(signal,luminance,reference,options)
% 每个像素的R/G/B共用缩放，保留代理通道比例；不逐通道截断。
    normalized=luminance/reference;
    target=options.peakLevel*(log1p(options.alpha*normalized)/log1p(options.alpha)).^options.gamma;
    gain=zeros(size(luminance)); active=luminance>0;
    gain(active)=target(active)./luminance(active);
    rgb=signal.*gain;
    compression=max(1,max(rgb,[],3)/options.peakLevel);
    rgb=rgb./compression;
end
