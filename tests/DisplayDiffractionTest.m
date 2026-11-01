classdef DisplayDiffractionTest < matlab.unittest.TestCase
    %DisplayDiffractionTest Analytical optics and image-formation regression tests.
    % Small rectangular pupils give known diffraction efficiencies and positions.
    % FFT bins and output bins coincide, so these tests do not hide window loss.
    methods (TestClassSetup)
        function addProjectSource(testCase)
            root = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'src')));
        end
    end

    methods (Test)
        function reflectionConservesKnownReflectance(testCase)
            cfg = smallConfig("reflection");
            cfg.featureIntensity = 0.36;
            cfg.backgroundIntensity = 0.36;

            result = runCOEReflection(cfg);

            testCase.verifyEqual(result.metrics.Throughput, 0.36, 'AbsTol', 2e-13);
            testCase.verifyEqual(sum(result.kernel(:)), 0.36, 'AbsTol', 2e-13);
            testCase.verifyLessThan(result.metrics.ParsevalAbsoluteError, 2e-13);
            testCase.verifyEqual(result.metrics.ZeroOrderRetention, 1, 'AbsTol', 2e-13);
        end

        function transmissionUsesIntensityNotItsSquare(testCase)
            cfg = smallConfig("udc");
            cfg.featureIntensity = 0.25;
            cfg.backgroundIntensity = 0.25;

            result = simulateDisplayDiffraction(cfg);

            testCase.verifyEqual(result.metrics.Throughput, 0.25, 'AbsTol', 2e-13);
            testCase.verifyEqual(sum(result.kernel(:)), 0.25, 'AbsTol', 2e-13);
            testCase.verifyEqual(sum(result.psf(:)), 1, 'AbsTol', 2e-13);
            testCase.verifyLessThan(result.metrics.ParsevalAbsoluteError, 2e-13);
        end

        function patternedCircularPupilHasAreaWeightedTransmission(testCase)
            cfg = smallConfig("udc");
            cfg.pupilRadius = 80e-6;
            cfg.featureIntensity = 0.65;
            cfg.backgroundIntensity = 0.02;

            result = simulateDisplayDiffraction(cfg);
            pupilIndices = result.pupil > 0;
            expected = mean(0.02 + 0.63*result.mask(pupilIndices));

            testCase.verifyEqual(result.metrics.Throughput, expected, 'AbsTol', 2e-13);
            testCase.verifyEqual(sum(result.kernel(:)), expected, 'AbsTol', 2e-13);
            testCase.verifyLessThan(result.metrics.ParsevalAbsoluteError, 2e-13);
        end

        function binaryHeightGivesSinglePassQuarterWave(testCase)
            cfg = smallConfig("udc");
            cfg.complexCoefficient = 1;
            cfg.heightMap = [zeros(64, 32), ones(64, 32)]*cfg.wavelengths/(4*cfg.deltaIndex);

            result = simulateDisplayDiffraction(cfg);

            % Half the pupil has phase pi/2: |(1+i)/2|^2 = 1/2.
            % Zero padding by two gives a clear-pupil centre bin of 1/4.
            testCase.verifyEqual(result.kernel(65, 65), 1/8, 'AbsTol', 2e-13);
            testCase.verifyEqual(result.metrics.Throughput, 1, 'AbsTol', 2e-13);
        end

        function binaryHeightCancelsDoublePassZeroOrder(testCase)
            cfg = smallConfig("reflection");
            cfg.complexCoefficient = 1;
            cfg.heightMap = [zeros(64, 32), ones(64, 32)]*cfg.wavelengths/(4*cfg.deltaIndex);

            result = runCOEReflection(cfg);

            % The same height now gives phase pi on half the pupil.
            testCase.verifyEqual(result.kernel(65, 65), 0, 'AbsTol', 2e-13);
            testCase.verifyEqual(result.baselineKernel(65, 65), 1/4, 'AbsTol', 2e-13);
            testCase.verifyEqual(result.metrics.Throughput, 1, 'AbsTol', 2e-13);
        end

        function uniformHeightOnlyChangesGlobalPhase(testCase)
            cfg = smallConfig("reflection");
            cfg.heightMap = 137e-9;

            result = runCOEReflection(cfg);

            testCase.verifyEqual(result.kernel, result.baselineKernel, 'AbsTol', 2e-13);
            testCase.verifyEqual(result.metrics.ZeroOrderRetention, 1, 'AbsTol', 2e-13);
        end

        function phaseGratingHasCorrectPositiveFirstOrder(testCase)
            cfg = smallConfig("reflection");
            cfg.complexCoefficient = repmat(exp(2i*pi*(0:63)/16), 64, 1);

            result = simulateDisplayDiffraction(cfg);
            [peak, location] = max(result.kernel(:));
            [row, col] = ind2sub(size(result.kernel), location);

            testCase.verifyEqual(row, 65);
            testCase.verifyEqual(col, 73);
            testCase.verifyEqual(result.coordinate(col), cfg.wavelengths/cfg.pitch, 'AbsTol', 2e-14);
            testCase.verifyEqual(peak, 1/4, 'AbsTol', 2e-13);
        end

        function wavelengthAndFocalLengthSetSensorOrderPosition(testCase)
            cfg = smallConfig("udc");
            cfg.complexCoefficient = repmat(exp(2i*pi*(0:63)/16), 64, 1);

            result = simulateDisplayDiffraction(cfg);
            [~, location] = max(result.kernel(:));
            [row, col] = ind2sub(size(result.kernel), location);

            testCase.verifyEqual(row, 65);
            testCase.verifyEqual(col, 73);
            testCase.verifyEqual(result.coordinate(col), ...
                cfg.wavelengths*cfg.focalLength/cfg.pitch, 'AbsTol', 2e-14);
        end

        function spectralWeightsAddIntensityIncoherently(testCase)
            cfg = smallConfig("reflection");
            cfg.wavelengths = [500, 500]*1e-9;
            cfg.spectralWeights = [1, 3];
            cfg.channelResponse = ones(2, 1);
            cfg.complexCoefficient = cat(3, ones(64), ...
                repmat(exp(2i*pi*(0:63)/16), 64, 1));

            result = simulateDisplayDiffraction(cfg);

            testCase.verifyEqual(result.kernel(65, 65), 1/16, 'AbsTol', 2e-13);
            testCase.verifyEqual(result.kernel(65, 73), 3/16, 'AbsTol', 2e-13);
            testCase.verifyEqual(result.channelSignal, 4*result.kernel, 'AbsTol', 2e-13);
            testCase.verifyEqual(sum(result.kernel(:)), 1, 'AbsTol', 2e-13);
        end

        function fixedExposurePreservesPerChannelLoss(testCase)
            cfg = smallConfig("udc");
            cfg.wavelengths = [500, 500, 500]*1e-9;
            cfg.spectralWeights = [1, 1, 1];
            cfg.channelResponse = eye(3);
            cfg.featureIntensity = [0.25, 0.36, 0.49];
            cfg.backgroundIntensity = cfg.featureIntensity;
            scene = repmat(reshape([0.2, 0.4, 0.6], 1, 1, 3), 9, 11, 1);

            [result, images] = runUDCImaging(cfg, scene, "linear");

            testCase.verifyEqual(result.channelThroughput, [0.25, 0.36, 0.49], 'AbsTol', 2e-13);
            testCase.verifyEqual(images.fixedExposure, ...
                scene.*reshape([0.25, 0.36, 0.49], 1, 1, 3), 'AbsTol', 2e-13);
            testCase.verifyEqual(images.brightnessMatched, scene, 'AbsTol', 2e-13);
            testCase.verifyEqual(images.clear, scene, 'AbsTol', 2e-13);
        end

        function rgbConvolutionMatchesDirectAsymmetricKernels(testCase)
            optics = asymmetricOptics();
            scene = reshape(mod(0:296, 31)/31, 9, 11, 3);
            expected = cat(3, directReplicate(scene(:, :, 1), optics.kernel(:, :, 1)), ...
                directReplicate(scene(:, :, 2), optics.kernel(:, :, 2)), ...
                directReplicate(scene(:, :, 3), optics.kernel(:, :, 3)));

            images = formUDCImage(scene, optics, "linear");

            testCase.verifyEqual(images.fixedExposure, expected, 'AbsTol', 2e-14);
            testCase.verifyEqual(images.clear, scene, 'AbsTol', 2e-14);
        end

        function hdrImpulsePreservesAmplitudeAndKernelOrientation(testCase)
            optics = asymmetricOptics();
            scene = zeros(9, 11, 3);
            scene(5, 6, :) = 1000;

            images = formUDCImage(scene, optics, "linear");

            testCase.verifyEqual(images.fixedExposure(4:6, 5:7, :), ...
                1000*optics.kernel, 'AbsTol', 2e-11);
            testCase.verifyEqual(max(images.inputLinear(:)), 1000);
            testCase.verifyGreaterThan(max(images.fixedExposure(:)), 1);
        end

        function srgbIsDecodedBeforeOpticalConvolution(testCase)
            optics = asymmetricOptics();
            scene = 0.5*ones(9, 11, 3);
            linearGray = ((0.5 + 0.055)/1.055)^2.4;
            totals = reshape(sum(sum(optics.kernel, 1), 2), 1, 1, 3);

            images = formUDCImage(scene, optics, "srgb");

            testCase.verifyEqual(images.inputLinear, linearGray*ones(9, 11, 3), 'AbsTol', 2e-14);
            testCase.verifyEqual(images.fixedExposure, ...
                linearGray*ones(9, 11, 3).*totals, 'AbsTol', 2e-14);
        end

        function zeroTransmissionHasNoNanOrSpuriousSignal(testCase)
            cfg = smallConfig("udc");
            cfg.featureIntensity = 0;
            cfg.backgroundIntensity = 0;

            [result, images] = runUDCImaging(cfg, ones(9, 11), "linear");

            testCase.verifyEqual(result.kernel, zeros(129));
            testCase.verifyEqual(result.psf, zeros(129));
            testCase.verifyEqual(result.mtf, zeros(129));
            testCase.verifyEqual(result.metrics.CapturedFraction, 0);
            testCase.verifyEqual(images.fixedExposure, zeros(9, 11));
            testCase.verifyEqual(images.brightnessMatched, zeros(9, 11));
            testCase.verifyEqual(images.clear, ones(9, 11), 'AbsTol', 2e-13);
        end

        function zeroReflectionHasUndefinedRetention(testCase)
            cfg = smallConfig("reflection");
            cfg.complexCoefficient = 0;

            result = runCOEReflection(cfg);

            testCase.verifyEqual(result.kernel, zeros(129));
            testCase.verifyTrue(isnan(result.metrics.ZeroOrderRetention));
        end

        function finiteWindowReportsTruncation(testCase)
            cfg = smallConfig("udc");
            cfg.complexCoefficient = 1;
            cfg.outputSize = 9;
            cfg.sensorPitch = 2e-6;
            cfg.paddingFactor = 8;

            testCase.verifyWarning(@() simulateDisplayDiffraction(cfg), ...
                'DisplayDiffraction:TruncatedEnergy');
        end

        function finiteWindowKernelRetainsPhysicalLoss(testCase)
            cfg = smallConfig("udc");
            cfg.complexCoefficient = 1;
            cfg.outputSize = 9;
            cfg.sensorPitch = 2e-6;
            cfg.paddingFactor = 8;
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture( ...
                'DisplayDiffraction:TruncatedEnergy'));

            result = simulateDisplayDiffraction(cfg);

            testCase.verifyEqual(result.metrics.Throughput, 1, 'AbsTol', 2e-13);
            testCase.verifyLessThan(result.metrics.CapturedThroughput, 0.98);
            testCase.verifyEqual(sum(result.kernel(:)), ...
                result.metrics.CapturedThroughput, 'AbsTol', 2e-13);
            testCase.verifyEqual(sum(result.psf(:)), 1, 'AbsTol', 2e-13);
        end

        function matchedSensorGridDoesNotWarnFromRoundoff(testCase)
            cfg = smallConfig("udc");

            testCase.verifyWarningFree(@() simulateDisplayDiffraction(cfg));
        end

        function normalizedMtfHasUnitDcAndBoundedMagnitude(testCase)
            cfg = smallConfig("udc");

            result = simulateDisplayDiffraction(cfg);

            testCase.verifyEqual(result.mtf(65, 65), 1, 'AbsTol', 2e-13);
            testCase.verifyLessThanOrEqual(max(result.mtf(:)), 1 + 2e-13);
            testCase.verifyEqual(result.spatialFrequency(65), 0, 'AbsTol', 2e-13);
        end

        function clippingPhysicalPupilIsRejected(testCase)
            cfg = smallConfig("udc");
            cfg.pupilRadius = 200e-6;

            testCase.verifyError(@() simulateDisplayDiffraction(cfg), 'DisplayDiffraction:PupilExtent');
        end

        function activeAmplitudeCoefficientIsRejected(testCase)
            cfg = smallConfig("udc");
            cfg.complexCoefficient = 1.01;

            testCase.verifyError(@() simulateDisplayDiffraction(cfg), 'DisplayDiffraction:PassiveCoefficient');
        end
    end
end

function cfg = smallConfig(mode)
    cfg = displayDiffractionDefaults(mode);
    cfg.periods = [4, 4];
    cfg.samplesPerPeriod = 16;
    cfg.paddingFactor = 2;
    cfg.pupilRadius = Inf; % Rectangular pupil for analytical checks.
    cfg.wavelengths = 500e-9;
    cfg.spectralWeights = 1;
    cfg.channelResponse = 1;
    cfg.outputSize = 129;
    cfg.sensorPitch = cfg.wavelengths*cfg.focalLength/(128*cfg.pitch/cfg.samplesPerPeriod);
    cfg.angleHalfRangeDeg = rad2deg(64*cfg.sensorPitch/cfg.focalLength);
end

function optics = asymmetricOptics()
    optics.config.mode = "udc";
    optics.kernel = cat(3, [0.01, 0.02, 0; 0.03, 0.12, 0.01; 0, 0.01, 0], ...
        [0, 0.02, 0.01; 0.03, 0.2, 0.04; 0.01, 0.02, 0.03], ...
        [0.03, 0, 0.02; 0.01, 0.25, 0.01; 0.02, 0.03, 0.03]);
    optics.psf = optics.kernel./sum(sum(optics.kernel, 1), 2);
    optics.clearKernel = repmat([0, 0, 0; 0, 1, 0; 0, 0, 0], 1, 1, 3);
end

function expected = directReplicate(image, kernel)
    % Independent direct spatial convolution, not the production FFT helper.
    padded = image([1, 1:end, end], [1, 1:end, end]);
    expected = conv2(padded, kernel, 'valid');
end
