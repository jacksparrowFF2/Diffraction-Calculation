classdef COESpotsTest < matlab.unittest.TestCase
    %COESpotsTest Shared ROIs, exact incident-energy normalization and exclusions.
    methods (TestClassSetup)
        function addModule(testCase)
            root = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'topic_coe_reflection', 'src')));
        end
    end
    methods (Test)
        function sharedROIsMeasureDifferentStatesWithoutNormalizingLoss(testCase)
            [a, b, spectrum, options] = syntheticSpots();
            spots = computeCOEDiffractionSpots(a, b, spectrum, options);
            rows = spots.metrics(spots.metrics.Group == "wavelength" & abs(spots.metrics.WavelengthNm - 500) < 1e-6, :);
            testCase.verifyEqual(height(rows), 4);
            left = rows(rows.CenterThetaXDeg < 0, :);
            right = rows(rows.CenterThetaXDeg > 0, :);
            testCase.verifyEqual(left.AbsoluteReflectancePercent, [20; 5], 'AbsTol', 1e-12);
            testCase.verifyEqual(right.AbsoluteReflectancePercent, [10; 30], 'AbsTol', 1e-12);
            testCase.verifyEqual(left.PeakThetaXDeg, [-1; -1], 'AbsTol', 1e-12);
            testCase.verifyEqual(right.PeakThetaXDeg, [1; 1], 'AbsTol', 1e-12);
            testCase.verifyEqual(rows.CenterThetaYDeg, zeros(4, 1), 'AbsTol', 1e-12);
            % 零级（20%）与边缘峰（5%）均不得进入统计。
            testCase.verifyGreaterThan(abs(rows.CenterThetaXDeg), 0.3);
            testCase.verifyLessThan(abs(rows.CenterThetaXDeg), 2);
            wide = spots.metrics(spots.metrics.Group == "radiant", :);
            testCase.verifyEqual(wide.AbsoluteReflectancePercent(wide.CenterThetaXDeg < 0), [35; 8.75], ...
                'AbsTol', 1e-12);
        end

        function minimumSeparationPreventsDoubleCounting(testCase)
            [a, ~, spectrum, options] = syntheticSpots();
            a.spectralEnergy(21, 32, :) = 0.09; % 离右侧亮斑 0.1° 的另一候选
            b = a;
            options.minSeparationDeg = 0;
            options.integrationHalfAngleDeg = 0.15;
            spots = computeCOEDiffractionSpots(a, b, spectrum, options);
            group = spots.groups(1);
            centers = group.centersDeg;
            for j = 1:size(centers, 1)
                others = [1:j-1, j+1:size(centers, 1)];
                testCase.verifyGreaterThan(hypot(centers(others, 1) - centers(j, 1), ...
                    centers(others, 2) - centers(j, 2)), 2*options.integrationHalfAngleDeg);
            end
            metrics = group.metrics;
            measured = sum(metrics.AbsoluteReflectancePercent(metrics.StateIndex == 1));
            testCase.verifyLessThanOrEqual(measured, 100*sum(a.spectralEnergy(:, :, 1), 'all') + 1e-12);
        end

        function zeroFieldReturnsEmptyTypedTable(testCase)
            [a, b, spectrum, options] = syntheticSpots();
            a.spectralEnergy(:) = 0; b.spectralEnergy(:) = 0;
            options.wavelengthNm = [];
            spots = computeCOEDiffractionSpots(a, b, spectrum, options);
            testCase.verifyEmpty(spots.metrics);
            testCase.verifyEqual(numel(spots.groups), 4);
            testCase.verifyTrue(ismember('AbsoluteReflectancePercent', spots.metrics.Properties.VariableNames));
        end

        function undersampledIntegrationIsRejected(testCase)
            [a, b, spectrum, options] = syntheticSpots();
            options.integrationHalfAngleDeg = 0.01;
            testCase.verifyError(@() computeCOEDiffractionSpots(a, b, spectrum, options), ...
                'COEReflection:SpotSampling');
        end
    end
end

function [a, b, spectrum, options] = syntheticSpots()
    a.config.wavelengths = [500, 600]*1e-9;
    a.coordinate = deg2rad(-2:0.1:2);
    a.spectralEnergy = zeros(41, 41, 2);
    a.spectralEnergy(21, 21, :) = 0.2;
    a.spectralEnergy(21, 11, :) = reshape([0.2, 0.4], 1, 1, 2);
    a.spectralEnergy(21, 31, :) = reshape([0.1, 0.2], 1, 1, 2);
    a.spectralEnergy(21, 41, :) = 0.05;
    b = a;
    b.spectralEnergy(21, 11, :) = reshape([0.05, 0.1], 1, 1, 2);
    b.spectralEnergy(21, 31, :) = reshape([0.3, 0.6], 1, 1, 2);
    spectrum = struct('energyWeights', [0.25; 0.75], 'photopicWeights', [0.6; 0.4]);
    options = struct('zeroOrderExclusionDeg', 0.3, 'integrationHalfAngleDeg', 0.15, ...
        'minSeparationDeg', 0.4, 'maxSpots', 8, 'minRelativePeak', 0.01, 'wavelengthNm', [500, 600]);
end
