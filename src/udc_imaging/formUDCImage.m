function images = formUDCImage(scene, optics, encoding)
%formUDCImage Linear intensity convolution; no clipping, tone mapping or noise.
% fixedExposure uses absolute kernels. brightnessMatched removes only overall
% captured loss. clear uses the same finite pupil without the display mask.
% Replicated scene edges suppress black padding; out-of-frame content is unknown.
    arguments
        scene {mustBeNumericOrLogical, mustBeNonempty}
        optics (1,1) struct
        encoding (1,1) string {mustBeMember(encoding, ["linear", "srgb"])} = "linear"
    end
    assert(optics.config.mode == "udc", 'DisplayDiffraction:Mode', 'UDC optics are required.');
    if isinteger(scene)
        assert(encoding == "srgb", 'DisplayDiffraction:Encoding', ...
            'Convert integer linear data to physical double values before simulation.');
        scene = double(scene)/double(intmax(class(scene)));
    else
        scene = double(scene);
    end
    validateattributes(scene, {'double'}, {'real', 'finite', 'nonnegative'});
    assert(ndims(scene) <= 3, 'DisplayDiffraction:SceneSize', 'Use a 2-D or channel image.');
    if encoding == "srgb"
        assert(all(scene(:) <= 1), 'DisplayDiffraction:Encoding', 'sRGB values must lie in [0,1].');
        low = scene <= 0.04045;
        scene(low) = scene(low)/12.92;
        scene(~low) = ((scene(~low) + 0.055)/1.055).^2.4;
    end
    channels = size(optics.kernel, 3);
    if size(scene, 3) == 1
        scene = repmat(scene, 1, 1, channels);
    end
    assert(size(scene, 3) == channels, 'DisplayDiffraction:SceneChannels', ...
        'Scene channel count must match channelResponse, or use a grayscale scene.');
    fixed = zeros(size(scene));
    matched = fixed;
    clearImage = fixed;
    for channel = 1:channels
        fixed(:, :, channel) = convolveReplicated(scene(:, :, channel), optics.kernel(:, :, channel));
        matched(:, :, channel) = convolveReplicated(scene(:, :, channel), optics.psf(:, :, channel));
        clearImage(:, :, channel) = convolveReplicated(scene(:, :, channel), optics.clearKernel(:, :, channel));
    end
    images.inputLinear = scene;
    images.clear = clearImage;
    images.fixedExposure = fixed;
    images.brightnessMatched = matched;
    images.boundary = "replicate";
end

function output = convolveReplicated(image, kernel)
    half = floor(size(kernel)/2);
    rows = min(max((1-half(1)):(size(image, 1)+half(1)), 1), size(image, 1));
    cols = min(max((1-half(2)):(size(image, 2)+half(2)), 1), size(image, 2));
    padded = image(rows, cols);
    % Reuse the project's existing linear FFT convolution; odd kernels avoid ambiguity.
    filtered = real(fftconv2_same(padded, kernel));
    output = max(0, filtered(half(1)+(1:size(image, 1)), half(2)+(1:size(image, 2))));
end
