function [withoutColor, withColor, mapping] = synthesizeCoeColor(withoutResult, withResult, spectrum, options)
%synthesizeCoeColor 按 D65 与 CIE XYZ 积分，将两种状态转换到共同曝光的 sRGB。
% 原始 XYZ/xy 不变；可选亮度映射只按 Y 计算统一比例缩放 XYZ，不逐通道增强。
    arguments
        withoutResult (1,1) struct
        withResult (1,1) struct
        spectrum (1,1) struct
        options (1,1) struct
    end
    validateattributes(options.exposureEV, {'numeric'}, {'scalar', 'real', 'finite'});
    validateattributes(options.fixedGain, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'});
    if ~isfield(options, 'toneMapping')
        options.toneMapping = "none";
    end
    if ~isfield(options, 'toneStrength')
        options.toneStrength = 3000;
    end
    if ~isfield(options, 'luminanceGamma')
        options.luminanceGamma = 0.75;
    end
    mustBeMember(options.toneMapping, ["none", "log"]);
    validateattributes(options.toneStrength, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
    validateattributes(options.luminanceGamma, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
    expected = spectrum.wavelengthNm(:)*1e-9;
    assert(numel(expected) == numel(withoutResult.config.wavelengths) ...
        && numel(expected) == numel(withResult.config.wavelengths) ...
        && all(abs(withoutResult.config.wavelengths(:) - expected) < 1e-15) ...
        && all(abs(withResult.config.wavelengths(:) - expected) < 1e-15), ...
        'COEReflection:ColorWavelengths', '色度数据与两种状态的波长网格必须一致。');
    withoutColor = integrateXYZ(withoutResult, spectrum);
    withColor = integrateXYZ(withResult, spectrum);
    % 转到每单位近轴立体角的 XYZ 密度，避免输出网格变密时图像自动变暗。
    reference = max([withoutColor.linearSRGBDensity(:); withColor.linearSRGBDensity(:); realmin]);
    gain = options.fixedGain;
    if gain == 0
        gain = 0.9/reference; % 两状态共同自动曝光；可用 fixedGain 固定跨方案曝光
    end
    gain = gain*2^options.exposureEV;
    mapping = options;
    mapping.gain = gain;
    mapping.referencePeak = reference;
    mapping.method = "Shared exposure and luminance-only XYZ scaling; peak compression for enhanced display";
    withoutColor = renderSRGB(withoutColor, gain, options);
    withColor = renderSRGB(withColor, gain, options);
end

function color = integrateXYZ(result, spectrum)
    gridSize = size(result.spectralEnergy, 1:2);
    signals = reshape(result.spectralEnergy, [], numel(spectrum.wavelengthNm));
    xyz = signals*spectrum.xyzWeights;
    color.XYZ = reshape(xyz, gridSize(1), gridSize(2), 3);
    total = sum(xyz, 2);
    xy = NaN(size(xyz, 1), 2);
    active = total > 0;
    xy(active, :) = xyz(active, 1:2)./total(active);
    color.xy = reshape(xy, gridSize(1), gridSize(2), 2);
    binSolidAngle = (result.coordinate(2) - result.coordinate(1))^2;
    color.XYZDensity = color.XYZ/binSolidAngle;
    xyzToRGB = [3.2404542, -1.5371385, -0.4985314; ...
        -0.9692660, 1.8760108, 0.0415560; 0.0556434, -0.2040259, 1.0572252];
    color.linearSRGBDensity = reshape((xyz/binSolidAngle)*xyzToRGB.', gridSize(1), gridSize(2), 3);
    color.binSolidAngle = binSolidAngle;
end

function color = renderSRGB(color, gain, options)
    rgb = color.linearSRGBDensity*gain;
    color.gamutClippedFraction = mean(any(rgb < -1e-7 | rgb > 1, 3), 'all');
    color.sRGB = encodeSRGB(rgb); % 保留原线性曝光版，便于与增强图比较
    scale = ones(size(rgb, 1:2));
    if options.toneMapping == "log"
        luminance = color.XYZDensity(:, :, 2)*gain;
        target = (log1p(options.toneStrength*luminance)/log1p(options.toneStrength)).^options.luminanceGamma;
        active = luminance > 0;
        scale(active) = target(active)./luminance(active);
        rgb = rgb.*scale;
        % 高光按整个 RGB 向量等比例缩小，不单独截断亮通道，以尽量保持色度。
        peakCompression = max(1, max(rgb, [], 3));
        rgb = rgb./peakCompression;
        scale = scale./peakCompression;
    end
    color.displayXYZ = color.XYZDensity.*(gain*scale);
    color.displayGamutClippedFraction = mean(any(rgb < -1e-7 | rgb > 1 + 1e-12, 3), 'all');
    color.displaySRGB = encodeSRGB(rgb);
end

function encoded = encodeSRGB(rgb)
    rgb = min(max(rgb, 0), 1); % 负线性 RGB 属于 sRGB 色域外，显示时仍需裁剪
    encoded = 12.92*rgb;
    nonlinear = rgb > 0.0031308;
    encoded(nonlinear) = 1.055*rgb(nonlinear).^(1/2.4) - 0.055;
end
