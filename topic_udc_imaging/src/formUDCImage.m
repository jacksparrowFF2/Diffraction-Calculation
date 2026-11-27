function images = formUDCImage(scene, optics, encoding, outputSize, frame)
%formUDCImage Linear intensity convolution; no clipping, tone mapping or noise.
% fixedExposure uses absolute kernels. brightnessMatched removes only overall
% captured loss. clear uses the same finite pupil without the display mask.
% Replicated scene edges suppress black padding; out-of-frame content is unknown.
    arguments
        scene {mustBeNumericOrLogical, mustBeNonempty}
        optics (1,1) struct
        encoding (1,1) string {mustBeMember(encoding, ["linear", "srgb"])} = "linear"
        outputSize double = []
        frame (1,1) struct = struct()
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
    if ~isempty(outputSize)
        validateattributes(outputSize,{'double'},{'row','numel',2,'integer','finite','positive'});
    end
    if ~isempty(outputSize) && ~isequal([size(scene,1),size(scene,2)],outputSize)
        scene = imresize(scene,outputSize,'bilinear','Antialiasing',true);
    end
    channels = size(optics.kernel, 3);
    if size(scene, 3) == 1
        scene = repmat(scene, 1, 1, channels);
    end
    assert(size(scene, 3) == channels, 'DisplayDiffraction:SceneChannels', ...
        'Scene channel count must match channelResponse, or use a grayscale scene.');
    if isfield(frame,'sizePixels')
        images = formFramedPreview(scene,optics,frame);
        return
    end
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

function images = formFramedPreview(scene,optics,frame)
    % Only the target and a finite PSF margin differ from the uniform background.
    background = frame.backgroundLinear;
    [rows,cols,channels] = size(scene);
    margin = floor([size(optics.kernel,1),size(optics.kernel,2)]/2);
    localSize = [rows,cols]+2*margin;
    previewSize = fliplr(frame.previewSizePixels);
    frameSize = fliplr(frame.sizePixels);
    rowWeights = frameOverlapWeights(localSize(1),frameSize(1),previewSize(1));
    colWeights = frameOverlapWeights(localSize(2),frameSize(2),previewSize(2)).';
    [input,clearImage,fixed,matched] = deal(zeros([previewSize,channels]));
    local = repmat(background,localSize);
    for channel = 1:channels
        local(margin(1)+(1:rows),margin(2)+(1:cols)) = scene(:,:,channel);
        clearKernel = optics.clearKernel(:,:,channel);
        fixedKernel = optics.kernel(:,:,channel);
        psf = optics.psf(:,:,channel);
        input(:,:,channel) = project(local,background);
        clearImage(:,:,channel) = project(convolveReplicated(local,clearKernel),background*sum(clearKernel,'all'));
        fixed(:,:,channel) = project(convolveReplicated(local,fixedKernel),background*sum(fixedKernel,'all'));
        matched(:,:,channel) = project(convolveReplicated(local,psf),background*sum(psf,'all'));
    end
    images = struct('inputLinear',input,'clear',clearImage,'fixedExposure',fixed, ...
        'brightnessMatched',matched,'boundary',"uniformBackground");

    function preview = project(image,base)
        preview = max(0,base+rowWeights*(image-base)*colWeights);
    end
end

function weights = frameOverlapWeights(localCount,frameCount,previewCount)
    % Average native detector pixels over the display pixels in linear light.
    firstOffset = floor((frameCount-localCount)/2);
    centers = firstOffset+(1:localCount)-(frameCount+1)/2;
    edges = linspace(-frameCount/2,frameCount/2,previewCount+1).';
    left = max(edges(1:end-1),centers-0.5);
    right = min(edges(2:end),centers+0.5);
    weights = sparse(max(0,right-left)/(frameCount/previewCount));
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
