classdef COESpectralTest < matlab.unittest.TestCase
    %COESpectralTest Physical normalization, angular integration and CIE colour checks.
    methods (TestClassSetup)
        function addModule(testCase)
            root = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'topic_coe_reflection', 'src')));
            addCoeAlgorithmFixture(testCase, root);
        end
    end
    methods (Test)
        function luminanceEnhancementPreservesXYZAndChromaticity(testCase)
            spectrum = loadCoeSpectrum(380:5:780);
            result.config.wavelengths = spectrum.wavelengthNm.'/1e9;
            result.coordinate = [-1, 0, 1]*1e-3;
            spatial = [0, 1e-7, 1e-5; 1e-4, 1e-3, 0.01; 0.1, 0.5, 1];
            spectralShape = 0.5 + exp(-((spectrum.wavelengthNm - 560)/70).^2);
            result.spectralEnergy = spatial.*reshape(spectralShape, 1, 1, []);
            options = struct('fixedGain', 0, 'exposureEV', 0, 'toneMapping', "log", ...
                'toneStrength', 3000, 'luminanceGamma', 0.75);
            [bright, identical] = synthesizeCoeColor(result, result, spectrum, options);
            options.toneMapping = "none";
            [original, ~] = synthesizeCoeColor(result, result, spectrum, options);
            testCase.verifyEqual(bright.XYZ, original.XYZ);
            testCase.verifyEqual(bright.xy, original.xy);
            testCase.verifyEqual(bright.linearSRGBDensity, original.linearSRGBDensity);
            testCase.verifyEqual(bright.sRGB, original.sRGB);
            testCase.verifyEqual(bright.displaySRGB, identical.displaySRGB);
            testCase.verifyGreaterThan(mean(bright.displaySRGB(2, 1, :)), mean(bright.sRGB(2, 1, :)));
            displayXYZ = reshape(bright.displayXYZ, [], 3);
            active = sum(displayXYZ, 2) > 0;
            xy = displayXYZ(active, 1:2)./sum(displayXYZ(active, :), 2);
            physicalXY = reshape(bright.xy, [], 2);
            testCase.verifyEqual(xy, physicalXY(active, :), 'AbsTol', 2e-13);
            % 本例色度位于 sRGB 色域内，编码后解码也应保持 RGB 分量比例。
            physicalRGB = reshape(bright.linearSRGBDensity, [], 3);
            decoded = reshape(decodeSRGB(bright.displaySRGB), [], 3);
            testCase.verifyGreaterThanOrEqual(min(physicalRGB(:)), 0);
            testCase.verifyEqual(decoded(active, :)./sum(decoded(active, :), 2), ...
                physicalRGB(active, :)./sum(physicalRGB(active, :), 2), 'AbsTol', 2e-13);
            testCase.verifyEqual(bright.displaySRGB(1, 1, :), zeros(1, 1, 3));
        end

        function officialD65HasExpectedWhitePointAndQuadrature(testCase)
            spectrum = loadCoeSpectrum(380:5:780);
            testCase.verifySize(spectrum.xyzWeights, [81, 3]);
            testCase.verifyEqual(sum(spectrum.quadratureNm), 400);
            testCase.verifyEqual(sum(spectrum.energyWeights), 1, 'AbsTol', 1e-14);
            testCase.verifyEqual(sum(spectrum.photopicWeights), 1, 'AbsTol', 1e-14);
            testCase.verifyEqual(spectrum.whiteXY, [0.3127, 0.3290], 'AbsTol', 1e-4);
            testCase.verifyEqual(spectrum.d65(spectrum.wavelengthNm == 560), 100);
        end

        function fresnelConstantAndMaterialMixConserveKnownReflectance(testCase)
            testCase.verifyEqual(normalIncidenceReflection(1, 1.5, 0), -0.2, 'AbsTol', 1e-14);
            testCase.verifyEqual(normalIncidenceReflection(1, 1, 0), 0);
            cfg = coeReflectionDefaults();
            cfg.periods = [4, 4]; cfg.samplesPerPeriod = 16; cfg.paddingFactor = 2;
            cfg.outputSize = 129; cfg.wavelengths = 550e-9; cfg.spectralWeights = 1; cfg.channelResponse = 1;
            cfg.mask = [zeros(64, 32), ones(64, 32)];
            cfg.materials = struct('enabled', true, 'incidentN', 1, 'featureN', 1.8, 'featureK', 0.05, ...
                'backgroundN', 1.5, 'backgroundK', 0);
            result = runCOEReflection(cfg);
            featureR = ((1 - 1.8)^2 + 0.05^2)/((1 + 1.8)^2 + 0.05^2);
            testCase.verifyEqual(result.metrics.Throughput, (featureR + 0.04)/2, 'AbsTol', 2e-13);
            cfg.complexCoefficient = 0.5;
            overridden = runCOEReflection(cfg);
            testCase.verifyEqual(overridden.metrics.Throughput, 0.25, 'AbsTol', 2e-13);
        end

        function angularIntegrationUsesIncidentEnergyAndAnnularRegion(testCase)
            result = toyResult();
            spectrum = struct('energyWeights', [0.25; 0.75], 'photopicWeights', [0.6; 0.4]);
            observation = struct('centerDeg', [0, 0], 'halfAngleDeg', 1, ...
                'innerHalfAngleDeg', 0, 'curveSamples', 11);
            angular = computeCoeAngularMetrics(result, observation, spectrum);
            testCase.verifyEqual(angular.spectralViewReflectance, [0.5; 0.5], 'AbsTol', 1e-14);
            testCase.verifyEqual(angular.radiantViewReflectance, 0.5, 'AbsTol', 1e-14);
            testCase.verifyEqual(angular.radiantTotalReflectance, 0.675, 'AbsTol', 1e-14);
            testCase.verifyEqual(angular.radiantCumulativeReflectance(1), 0);
            testCase.verifyEqual(angular.radiantCumulativeReflectance(end), 0.5, 'AbsTol', 1e-14);
            testCase.verifyGreaterThanOrEqual(diff(angular.radiantCumulativeReflectance), 0);
            observation.innerHalfAngleDeg = 0.1;
            ring = computeCoeAngularMetrics(result, observation, spectrum);
            testCase.verifyEqual(ring.radiantViewReflectance, 0.15, 'AbsTol', 1e-14);
            testCase.verifyEqual(ring.photopicViewReflectance, 0.22, 'AbsTol', 1e-14);
            observation.centerDeg = [0.5, 0]; observation.halfAngleDeg = 0.5;
            shifted = computeCoeAngularMetrics(result, observation, spectrum);
            testCase.verifyEqual(shifted.radiantViewReflectance, 0.5, 'AbsTol', 1e-14);
            observation.halfAngleDeg = 1;
            testCase.verifyError(@() computeCoeAngularMetrics(result, observation, spectrum), ...
                'COEReflection:ObservationExtent');
        end

        function neutralReflectorRetainsXYZAndCommonExposure(testCase)
            spectrum = loadCoeSpectrum(380:5:780);
            result.config.wavelengths = spectrum.wavelengthNm.'/1e9;
            result.coordinate = [-1, 0, 1]*1e-3;
            result.spectralEnergy = ones(3, 3, 81)*(0.25/9);
            darker = result;
            darker.spectralEnergy = result.spectralEnergy/2;
            options = struct('fixedGain', 1e-6, 'exposureEV', 0);
            [bright, dim] = synthesizeCoeColor(result, darker, spectrum, options);
            testCase.verifyEqual(sum(bright.XYZ(:, :, 2), 'all'), 0.25, 'AbsTol', 1e-13);
            testCase.verifyEqual(dim.XYZ, bright.XYZ/2, 'AbsTol', 1e-13);
            testCase.verifyEqual(bright.xy(2, 2, :), reshape(spectrum.whiteXY, 1, 1, 2), 'AbsTol', 1e-13);
            decoded = decodeSRGB(bright.sRGB);
            testCase.verifyEqual(decodeSRGB(dim.sRGB), decoded/2, 'AbsTol', 1e-13);
            testCase.verifyEqual(bright.gamutClippedFraction, 0);
            darker.spectralEnergy(:) = 0;
            [~, black] = synthesizeCoeColor(result, darker, spectrum, options);
            testCase.verifyEqual(black.sRGB, zeros(3, 3, 3));
            testCase.verifyTrue(all(isnan(black.xy(:))));
        end
    end
end

function result = toyResult()
    result.config.wavelengths = [500, 600]*1e-9;
    result.coordinate = deg2rad([-1, 0, 1]);
    result.spectralEnergy = zeros(3, 3, 2);
    result.spectralEnergy(2, 2, :) = reshape([0.2, 0.4], 1, 1, 2);
    result.spectralEnergy(2, 3, :) = reshape([0.3, 0.1], 1, 1, 2);
    result.spectralEnergy(1, 1, :) = reshape([0.1, 0.2], 1, 1, 2);
    result.metrics = table([0.6; 0.7], [0.6; 0.7], ...
        'VariableNames', {'Throughput', 'CapturedThroughput'});
end

function linear = decodeSRGB(encoded)
    linear = encoded/12.92;
    high = encoded > 0.04045;
    linear(high) = ((encoded(high) + 0.055)/1.055).^2.4;
end
