function result = simulateDisplayDiffraction(cfg)
%simulateDisplayDiffraction Scalar, normal-incidence reflection or pupil-plane UDC.
% Reflection uses the aperture module's chirp-convolution Fourier quadrature,
% including exp(i*k*z)/(i*lambda*z), without receiver quadratic phase.
% UDC retains its sampled FFT and energy-bin resampling.
% Wavelengths add INCOHERENTLY after squaring; r/t remain complex amplitudes.
% Output kernels are integrated energy per output bin, relative to an ideal
% clear pupil. psf is unit-sum shape; kernel retains loss AND finite-window loss.
    arguments
        cfg (1,1) struct
    end
    validateConfiguration(cfg);
    [mask, x, y] = buildGeometry(cfg);
    [xx, yy] = meshgrid(x, y);
    pupil = double(xx.^2 + yy.^2 <= cfg.pupilRadius^2);
    incidentEnergy = sum(pupil(:).^2);
    assert(incidentEnergy > 0, 'DisplayDiffraction:EmptyPupil', ...
        'The sampled pupil is empty. Increase pupilRadius or sampling density.');
    gridSize = size(mask);
    dx = cfg.pitch/cfg.samplesPerPeriod;
    assert(isinf(cfg.pupilRadius) || 2*cfg.pupilRadius <= min(gridSize)*dx, ...
        'DisplayDiffraction:PupilExtent', ...
        'The geometry must cover the complete pupil diameter. Increase periods or reduce pupilRadius.');
    height = expandMap(cfg.heightMap, gridSize, 'heightMap');
    validateattributes(height, {'numeric'}, {'real', 'finite'});
    fftSize = 2^nextpow2(max(gridSize)*cfg.paddingFactor);
    if cfg.mode == "reflection"
        fftSize = max(gridSize); % Only the independent Parseval FFT; no output padding.
    end
    frequency = (-fftSize/2:fftSize/2-1)/(fftSize*dx);
    outputCoordinate = outputGrid(cfg);
    binStep = outputCoordinate(2) - outputCoordinate(1);
    nLambda = numel(cfg.wavelengths);
    energy = zeros(cfg.outputSize, cfg.outputSize, nLambda);
    clearEnergy = energy;
    throughput = zeros(1, nLambda);
    captured = throughput;
    clearCaptured = throughput;
    parsevalError = throughput;
    zeroOrder = NaN(1, nLambda);
    outsideAngle = zeroOrder;
    nativeStep = throughput;
    receiverSamples = zeros(1, nLambda);
    receiverIrradiance = energy;
    phasePasses = 1 + double(cfg.mode == "reflection");
    [outX, outY] = meshgrid(outputCoordinate);
    metricMask = hypot(outX, outY) <= deg2rad(cfg.metricAngleDeg);

    for k = 1:nLambda
        lambda = cfg.wavelengths(k);
        coefficient = opticalCoefficient(cfg, mask, k);
        phase = phasePasses*2*pi*cfg.deltaIndex*height/lambda;
        field = pupil.*coefficient.*exp(1i*phase);
        if cfg.mode == "reflection"
            % Unpadded FFT is only an independent discrete Parseval diagnostic.
            nativeEnergy = transformEnergy(field, max(gridSize), incidentEnergy);
        else
            nativeEnergy = transformEnergy(field, fftSize, incidentEnergy);
        end
        throughput(k) = sum(abs(field(:)).^2)/incidentEnergy;
        parsevalError(k) = abs(sum(nativeEnergy(:)) - throughput(k));
        nativeCoordinate = lambda*frequency;
        if cfg.mode == "udc"
            nativeCoordinate = cfg.focalLength*nativeCoordinate;
        end
        nativeStep(k) = nativeCoordinate(2) - nativeCoordinate(1);
        if cfg.mode == "reflection"
            incidentArea = incidentEnergy*dx^2;
            [energy(:, :, k), intensityRatio, receiverSamples(k)] = ...
                propagateCoeReceiver(field, x, y, outputCoordinate, lambda, cfg, incidentArea);
            receiverIrradiance(:, :, k) = cfg.incidentIrradiance*intensityRatio;
            clearEnergy(:, :, k) = propagateCoeReceiver(pupil, x, y, ...
                outputCoordinate, lambda, cfg, incidentArea);
        else
            weights = binOverlap(nativeCoordinate, outputCoordinate, binStep);
            energy(:, :, k) = weights*nativeEnergy*weights.';
            clearNative = transformEnergy(pupil, fftSize, incidentEnergy);
            clearEnergy(:, :, k) = weights*clearNative*weights.';
        end
        captured(k) = sum(energy(:, :, k), 'all');
        clearCaptured(k) = sum(clearEnergy(:, :, k), 'all');
        if cfg.mode == "reflection" && throughput(k) > 0
            zeroRadius = cfg.zeroOrderRadiusFactor*lambda/cfg.pitch;
            zeroMask = hypot(outX, outY) <= zeroRadius;
            zeroOrder(k) = sum(energy(:, :, k).*zeroMask, 'all')/throughput(k);
            inside = sum(energy(:, :, k).*metricMask, 'all');
            outsideAngle(k) = max(0, 1 - inside/throughput(k));
        end
    end
    response = cfg.channelResponse.*cfg.spectralWeights(:);
    channelWeight = sum(response, 1);
    channelWeight(channelWeight == 0) = 1;
    response = response./channelWeight;
    channels = size(response, 2);
    kernel = reshape(reshape(energy, [], nLambda)*response, ...
        cfg.outputSize, cfg.outputSize, channels);
    clearKernel = reshape(reshape(clearEnergy, [], nLambda)*response, ...
        cfg.outputSize, cfg.outputSize, channels);
    channelCaptured = reshape(sum(sum(kernel, 1), 2), 1, []);
    channelThroughput = throughput*response;
    psf = normalizeKernels(kernel);
    clearPsf = normalizeKernels(clearKernel);
    result.config = cfg;
    result.mask = mask;
    result.pupil = pupil;
    result.x = x;
    result.y = y;
    result.coordinate = outputCoordinate;
    result.coordinateUnit = "radian (paraxial)";
    if cfg.mode == "udc"
        result.coordinateUnit = "m (sensor)";
    end
    result.phasePasses = phasePasses;
    result.spectralEnergy = energy;
    result.kernel = kernel;
    result.channelSignal = kernel.*reshape(sum(cfg.channelResponse.*cfg.spectralWeights(:), 1), 1, 1, []);
    result.psf = psf;
    result.clearKernel = clearKernel;
    result.clearPsf = clearPsf;
    result.metrics = table(cfg.wavelengths(:)*1e9, throughput(:), captured(:), ...
        safeRatio(captured(:), throughput(:)), clearCaptured(:), parsevalError(:), ...
        zeroOrder(:), outsideAngle(:), 'VariableNames', {'WavelengthNm', ...
        'Throughput', 'CapturedThroughput', 'CapturedFraction', 'ClearCapturedFraction', ...
        'ParsevalAbsoluteError', 'ZeroOrderFraction', 'OutsideAngleFraction'});
    result.channelThroughput = channelThroughput;
    result.channelCapturedThroughput = channelCaptured;
    result.nativeStep = nativeStep;
    result.fftSize = fftSize;
    if cfg.mode == "reflection"
        result.propagationMethod = "Fraunhofer chirp-convolution Fourier integral";
        result.receiverCoordinate = cfg.propagationDistance*outputCoordinate; % m, X=z*theta.
        result.receiverPixelArea = (cfg.propagationDistance*binStep)^2;
        result.incidentArea = incidentEnergy*dx^2;
        result.incidentPower = cfg.incidentIrradiance*result.incidentArea;
        result.spectralIrradiance = receiverIrradiance; % Pixel-averaged W/m^2.
        result.receiverSamplesPerBin = receiverSamples;
        result.fraunhoferNumber = (max(gridSize)*dx)^2 ...
            /(min(cfg.wavelengths)*cfg.propagationDistance);
        if result.fraunhoferNumber > 0.5
            warning('COEReflection:FraunhoferDistance', ...
                'z=%.4g m 小于常用远场尺度 2*D^2/lambda_min=%.4g m；本计算仍使用夫琅禾费近似，请增加距离或使用近场传播模型。', ...
                cfg.propagationDistance, 2*(max(gridSize)*dx)^2/min(cfg.wavelengths));
        end
    end
    if cfg.mode == "udc"
        [result.mtf, result.spatialFrequency] = computeMtf(psf, cfg.sensorPitch);
    end
    capturedFraction = result.metrics.CapturedFraction;
    truncated = find(capturedFraction < 0.98 & throughput(:) > 0);
    if ~isempty(truncated)
        % 只列出未达到阈值的波长；两种能量分母明确区分，不修改物理数组。
        message = sprintf('以下 %d 个波长在输出窗口内捕获的反射能量不足 98%%。\n', numel(truncated));
        if cfg.mode == "reflection"
            message = [message, sprintf('当前观察角域：θx、θy 均为 [-%.3f°, +%.3f°]。\n', ...
                cfg.angleHalfRangeDeg, cfg.angleHalfRangeDeg)];
        end
        message = [message, sprintf('波长/nm\t窗口内/总反射\t窗口外/总反射\t窗口内/入射\n')];
        for index = truncated(:).'
            message = [message, sprintf('%.3f\t%.3f%%\t\t%.3f%%\t\t%.3f%%\n', ...
                result.metrics.WavelengthNm(index), 100*capturedFraction(index), ...
                100*(1-capturedFraction(index)), 100*captured(index))]; %#ok<AGROW>
        end
        message = [message, '这是观察窗口截断提示，不代表能量不守恒；可查看 result.metrics.CapturedFraction。', ...
            '如需覆盖更完整的远场，可增大 cfg.angleHalfRangeDeg，并相应提高 cfg.outputSize 保持角度分辨率。'];
        warning('DisplayDiffraction:TruncatedEnergy', '%s', message);
    end
    % Computing spacing by subtraction can round a matched grid slightly up.
    if cfg.mode == "udc" && any(nativeStep > cfg.sensorPitch*(1 + 1e-12))
        warning('DisplayDiffraction:CoarsePSF', ...
            'Native PSF spacing exceeds sensorPitch. Increase paddingFactor and check convergence.');
    end
end

function validateConfiguration(cfg)
    mustBeMember(cfg.mode, ["reflection", "udc"]);
    names = {'pitch', 'radius', 'aspectRatio', 'focalLength', 'sensorPitch', 'angleHalfRangeDeg'};
    for k = 1:numel(names)
        validateattributes(cfg.(names{k}), {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
    end
    validateattributes(cfg.periods, {'numeric'}, {'row', 'numel', 2, 'integer', 'positive'});
    assert(all(mod(cfg.periods, 2) == 0), 'DisplayDiffraction:Periods', ...
        'periods must contain even counts for the existing 2-by-2 cell.');
    validateattributes(cfg.samplesPerPeriod, {'numeric'}, {'scalar', 'integer', '>=', 4});
    validateattributes(cfg.paddingFactor, {'numeric'}, {'scalar', 'real', 'finite', '>=', 1});
    validateattributes(cfg.outputSize, {'numeric'}, {'scalar', 'integer', '>=', 3});
    assert(mod(cfg.outputSize, 2) == 1, 'DisplayDiffraction:OutputSize', 'outputSize must be odd.');
    validateattributes(cfg.pupilRadius, {'numeric'}, {'scalar', 'real', 'positive', 'nonnan'});
    validateattributes(cfg.deltaIndex, {'numeric'}, {'scalar', 'real', 'finite'});
    validateattributes(cfg.zeroOrderRadiusFactor, {'numeric'}, {'scalar', 'real', '>', 0, '<', 0.5});
    validateattributes(cfg.metricAngleDeg, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
    validateattributes(cfg.cellAnglesDeg, {'numeric'}, {'size', [2, 2], 'real', 'finite'});
    validateattributes(cfg.wavelengths, {'numeric'}, {'row', 'real', 'finite', 'positive', 'nonempty'});
    nLambda = numel(cfg.wavelengths);
    validateattributes(cfg.spectralWeights, {'numeric'}, {'vector', 'numel', nLambda, 'real', 'finite', 'nonnegative'});
    validateattributes(cfg.channelResponse, {'numeric'}, {'2d', 'real', 'finite', 'nonnegative', 'nonempty'});
    assert(size(cfg.channelResponse, 1) == nLambda, 'DisplayDiffraction:Spectrum', ...
        'channelResponse needs one row per wavelength.');
    assert(any(cfg.channelResponse.*cfg.spectralWeights(:) > 0, 'all'), ...
        'DisplayDiffraction:Spectrum', 'At least one channel needs nonzero spectral weight.');
    if cfg.mode == "reflection"
        validateattributes(cfg.propagationDistance, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
        validateattributes(cfg.incidentIrradiance, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'});
        validateattributes(cfg.receiverBinSamples, {'numeric'}, {'scalar', 'integer', 'positive', 'odd'});
        dx = cfg.pitch/cfg.samplesPerPeriod;
        assert(deg2rad(cfg.angleHalfRangeDeg)*(1 + 1/(cfg.outputSize-1)) ...
            <= min(cfg.wavelengths)/(2*dx)*(1 + 1e-12), 'COEReflection:ReceiverAliasing', ...
            '输出像素边界超过物面采样的 Nyquist 角域；请提高 samplesPerPeriod 或缩小 angleHalfRangeDeg。');
        assert(cfg.angleHalfRangeDeg <= 10, 'DisplayDiffraction:ParaxialRange', ...
            'This paraxial model supports angleHalfRangeDeg <= 10 degrees.');
        assert(cfg.metricAngleDeg <= cfg.angleHalfRangeDeg, 'DisplayDiffraction:MetricAngle', ...
            'metricAngleDeg must lie inside the output window.');
    end
end

function [mask, x, y] = buildGeometry(cfg)
    n = cfg.samplesPerPeriod;
    dx = cfg.pitch/n;
    local = ((0:n-1) - floor(n/2))*dx; % Excludes duplicate pitch endpoint.
    x = ((0:cfg.periods(1)*n-1) - floor(cfg.periods(1)*n/2))*dx;
    y = ((0:cfg.periods(2)*n-1) - floor(cfg.periods(2)*n/2))*dx;
    if isempty(cfg.mask)
        [cx, cy] = meshgrid(local);
        unit = zeros(2*n);
        for row = 1:2
            for col = 1:2
                unit((row-1)*n+(1:n), (col-1)*n+(1:n)) = customPattern( ...
                    cfg.cellAnglesDeg(row, col), cx, cy, cfg.radius*cfg.aspectRatio, cfg.radius, 0, 0);
            end
        end
        mask = repmat(unit, cfg.periods(2)/2, cfg.periods(1)/2);
    else
        mask = cfg.mask;
        validateattributes(mask, {'numeric', 'logical'}, ...
            {'real', 'finite', 'size', [numel(y), numel(x)], '>=', 0, '<=', 1});
        mask = double(mask);
    end
end

function coefficient = opticalCoefficient(cfg, mask, k)
    if isempty(cfg.complexCoefficient)
        if isfield(cfg, 'coeStack') && cfg.coeStack.enabled
            assert(isequal(size(cfg.coeStack.regionMap), size(mask)), ...
                'COEReflection:RegionSize', '膜堆区域图必须与物面采样相同。');
            values = cfg.coeStack.coefficientLookup(:, k);
            coefficient = reshape(values(double(cfg.coeStack.regionMap(:)) + 1), size(mask));
        elseif isfield(cfg, 'materials') && cfg.materials.enabled
            % 两种材料的常数 n/k：法向单界面 Fresnel 复系数，包含材料反射相位。
            feature = normalIncidenceReflection(cfg.materials.incidentN, ...
                cfg.materials.featureN, cfg.materials.featureK);
            background = normalIncidenceReflection(cfg.materials.incidentN, ...
                cfg.materials.backgroundN, cfg.materials.backgroundK);
            assert(isscalar(feature) && isscalar(background), 'COEReflection:ConstantNK', ...
                '当前材料模型要求每种材料的 n、k 为常数标量。');
            coefficient = background + (feature - background)*mask;
        else
            feature = spectralValue(cfg.featureIntensity, k, numel(cfg.wavelengths));
            background = spectralValue(cfg.backgroundIntensity, k, numel(cfg.wavelengths));
            % mask interpolates material area only; prefer binary mask or complexCoefficient.
            coefficient = sqrt(background) + (sqrt(feature) - sqrt(background))*mask;
        end
    else
        value = cfg.complexCoefficient;
        assert(ndims(value) <= 3 && (size(value, 3) == 1 || ...
            size(value, 3) == numel(cfg.wavelengths)), 'DisplayDiffraction:CoefficientSize', ...
            'complexCoefficient must be scalar, Ny-by-Nx or Ny-by-Nx-by-wavelength.');
        coefficient = expandMap(value(:, :, min(k, size(value, 3))), size(mask), 'complexCoefficient');
        validateattributes(coefficient, {'numeric'}, {'finite'});
        assert(all(abs(coefficient(:)) <= 1 + 1e-12), 'DisplayDiffraction:PassiveCoefficient', ...
            'Use a passive power-normalized amplitude coefficient with abs(r/t) <= 1.');
    end
end

function value = spectralValue(values, k, count)
    validateattributes(values, {'numeric'}, {'vector', 'real', 'finite', '>=', 0, '<=', 1});
    assert(isscalar(values) || numel(values) == count, 'DisplayDiffraction:Spectrum', ...
        'Intensity coefficients must be scalar or one value per wavelength.');
    value = values(min(k, numel(values)));
end

function map = expandMap(value, gridSize, name)
    assert(isnumeric(value) && (isscalar(value) || isequal(size(value), gridSize)), ...
        'DisplayDiffraction:MapSize', '%s must be scalar or match the sampled geometry.', name);
    map = value.*ones(gridSize);
end

function energy = transformEnergy(field, fftSize, incidentEnergy)
    % Placement only introduces a Fourier phase ramp, which cancels in intensity.
    spectrum = fftshift(fft2(field, fftSize, fftSize));
    energy = abs(spectrum).^2/(fftSize^2*incidentEnergy);
end

function coordinate = outputGrid(cfg)
    half = (cfg.outputSize - 1)/2;
    if cfg.mode == "udc"
        coordinate = (-half:half)*cfg.sensorPitch;
    else
        coordinate = linspace(-deg2rad(cfg.angleHalfRangeDeg), ...
            deg2rad(cfg.angleHalfRangeDeg), cfg.outputSize);
    end
end

function weights = binOverlap(native, target, targetStep)
    % Conserves bin energy under resampling, assuming constant density per FFT bin.
    nativeStep = native(2) - native(1);
    left = max(target(:) - targetStep/2, native - nativeStep/2);
    right = min(target(:) + targetStep/2, native + nativeStep/2);
    weights = sparse(max(0, right - left)/nativeStep);
end

function normalized = normalizeKernels(kernel)
    total = sum(sum(kernel, 1), 2);
    total(total == 0) = 1;
    normalized = kernel./total;
end

function ratio = safeRatio(numerator, denominator)
    ratio = zeros(size(numerator));
    active = denominator > 0;
    ratio(active) = numerator(active)./denominator(active);
end

function [mtf, frequency] = computeMtf(psf, sensorPitch)
    mtf = zeros(size(psf));
    for channel = 1:size(psf, 3)
        mtf(:, :, channel) = abs(fftshift(fft2(ifftshift(psf(:, :, channel)))));
    end
    n = size(psf, 1);
    frequency = (-floor(n/2):floor(n/2))/(n*sensorPitch*1e3); % cycles/mm
end
