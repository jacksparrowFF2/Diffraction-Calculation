classdef COEMultilayerTest < matlab.unittest.TestCase
    %COEMultilayerTest Analytic film limits, absorber footprint and layer input contracts.
    methods (TestClassSetup)
        function addModule(testCase)
            root = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'topic_coe_reflection', 'src')));
        end
    end
    methods (Test)
        function bareMirrorAbsorberAndCommonReferencePlane(testCase)
            wavelengths = [450, 550, 650]*1e-9;
            mirror = reflectCOEFilmStack(wavelengths, 1, zeros(0, 3), zeros(1, 0), -1, 0);
            absorber = reflectCOEFilmStack(wavelengths, 1, zeros(0, 3), zeros(1, 0), 0, 0);
            testCase.verifyEqual(mirror, -ones(1, 3));
            testCase.verifyEqual(absorber, zeros(1, 3));
            airLayer = reflectCOEFilmStack(wavelengths, 1, ones(1, 3), 100e-9, -1, 500e-9);
            bareAtReference = reflectCOEFilmStack(wavelengths, 1, zeros(0, 3), zeros(1, 0), -1, 500e-9);
            testCase.verifyEqual(airLayer, bareAtReference, 'AbsTol', 2e-13);
        end

        function quarterWaveFilmCancelsKnownInterface(testCase)
            lambda = 550e-9;
            filmN = sqrt(1.5);
            depth = lambda/(4*filmN);
            coefficient = reflectCOEFilmStack(lambda, 1, filmN, depth, -0.2, depth);
            testCase.verifyEqual(coefficient, 0, 'AbsTol', 2e-13);
            zeroDepth = reflectCOEFilmStack(lambda, 1, filmN, 0, -0.2, 0);
            testCase.verifyEqual(zeroDepth, -0.2, 'AbsTol', 2e-13);
        end

        function losslessStackOnMirrorConservesUnitReflectance(testCase)
            wavelengths = (380:5:780)*1e-9;
            films = [1.7; 2.1].*ones(2, numel(wavelengths));
            depths = [120, 250]*1e-9;
            coefficient = reflectCOEFilmStack(wavelengths, 1, films, depths, -1, sum(depths));
            testCase.verifyEqual(abs(coefficient).^2, ones(size(wavelengths)), 'AbsTol', 2e-13);
            lossy = reflectCOEFilmStack(wavelengths, 1, 1.7 + 0.2i, 600e-9, -1, 600e-9);
            testCase.verifyLessThan(abs(lossy).^2, 1);
        end

        function opaqueFilmHidesUnderlyingTermination(testCase)
            wavelengths = [450, 550, 650]*1e-9;
            film = 1.8 + 5i;
            mirror = reflectCOEFilmStack(wavelengths, 1, film, 1e-6, -1, 1e-6);
            absorber = reflectCOEFilmStack(wavelengths, 1, film, 1e-6, 0, 1e-6);
            expected = (1 - film)/(1 + film);
            testCase.verifyEqual(mirror, expected*ones(1, 3), 'AbsTol', 2e-13);
            testCase.verifyEqual(absorber, mirror, 'AbsTol', 2e-13);
        end

        function absorbingTwoLayerStackMatchesCharacteristicMatrix(testCase)
            wavelengths = [450, 550, 650]*1e-9;
            indices = [1.7 + 0.05i; 2.0 + 0.1i];
            depths = [120, 250]*1e-9;
            expected = zeros(size(wavelengths));
            for k = 1:numel(wavelengths)
                matrix = eye(2);
                for layer = 1:2
                    delta = 2*pi*indices(layer)*depths(layer)/wavelengths(k);
                    transfer = [cos(delta), -1i*sin(delta)/indices(layer); ...
                        -1i*indices(layer)*sin(delta), cos(delta)];
                    matrix = matrix*transfer;
                end
                admittance = (matrix(2, 1) + matrix(2, 2))/(matrix(1, 1) + matrix(1, 2));
                expected(k) = (1 - admittance)/(1 + admittance);
            end
            actual = reflectCOEFilmStack(wavelengths, 1, indices, depths, 0, sum(depths));
            testCase.verifyEqual(actual, expected, 'AbsTol', 2e-13);
        end

        function expandedBMDoesNotExpandAnodeReflectance(testCase)
            cfg = smallConfig(); stack = coeStackDefaults();
            for name = ["BM", "R", "G", "B"]
                stack.materials.(name).thickness = 0;
            end
            [model, geometry] = prepareCOEStack(cfg, stack);
            testCase.verifyGreaterThan(nnz(geometry.expansionAbsorberMask), 0);
            testCase.verifyFalse(any(geometry.anodeMask & geometry.expansionAbsorberMask, 'all'));
            testCase.verifyFalse(any(double(geometry.cfRMask) + double(geometry.cfGMask) + double(geometry.cfBMask) > 1, 'all'));
            for k = 1:3
                lookup = model.coefficientLookup(:, k);
                coefficient = reshape(lookup(double(model.regionMap(:)) + 1), size(model.regionMap));
                testCase.verifyEqual(abs(coefficient), double(geometry.anodeMask), 'AbsTol', 2e-13);
            end
            cfg.coeStack = model; cfg.mask = double(geometry.anodeMask);
            result = runCOEReflection(cfg);
            testCase.verifyEqual(result.metrics.Throughput, mean(geometry.anodeMask, 'all')*ones(3, 1), ...
                'AbsTol', 2e-13);
            testCase.verifyEqual(result.kernel, result.baselineKernel, 'AbsTol', 2e-13);
        end

        function layerMATAndImagesRoundTripAndCSVNK(testCase)
            temporary = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture());
            cfg = smallConfig(); stack = coeStackDefaults();
            [expectedModel, expected] = prepareCOEStack(cfg, stack);
            fields = ["anodeMask", "bmOpenMask", "cfRMask", "cfGMask", "cfBMask"];
            data = struct();
            for name = fields
                data.(name) = expected.(name);
            end
            save(fullfile(temporary.Folder, 'coe_layers.mat'), '-struct', 'data');
            stack.geometry.source = "mat";
            [fromMAT, actual] = prepareCOEStack(cfg, stack, temporary.Folder);
            testCase.verifyEqual(actual.anodeMask, expected.anodeMask);
            testCase.verifyEqual(fromMAT.coefficientLookup, expectedModel.coefficientLookup);
            stack.geometry.source = "images";
            imwrite(data.anodeMask, fullfile(temporary.Folder, 'anode.png'));
            imwrite(~data.bmOpenMask, fullfile(temporary.Folder, 'bm.png'));
            for j = 1:3
                imwrite(data.(fields(j+2)), fullfile(temporary.Folder, stack.geometry.cfFiles(j)));
            end
            stack.geometry.bmFile = "bm.png"; stack.geometry.bmImageIsOpening = false;
            writematrix([450, 1.9, 0.35; 550, 1.8, 0.35; 650, 1.7, 0.35], fullfile(temporary.Folder, 'bm_nk.csv'));
            stack.materials.BM.nkFile = "bm_nk.csv";
            [fromImages, actual] = prepareCOEStack(cfg, stack, temporary.Folder);
            testCase.verifyEqual(actual.cfLabel, expected.cfLabel);
            testCase.verifyEqual(fromImages.materialTable.n(1:3), [1.9; 1.8; 1.7]);
        end

        function CFOverlapAndMissingCoverageAreRejected(testCase)
            cfg = smallConfig(); stack = coeStackDefaults();
            stack.demo.cfExpansion = 20e-6;
            testCase.verifyError(@() prepareCOEStack(cfg, stack), 'COEReflection:CFOverlap');
            stack = coeStackDefaults();
            [~, data] = prepareCOEStack(cfg, stack);
            temporary = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture());
            data.cfGMask(:) = false;
            save(fullfile(temporary.Folder, 'coe_layers.mat'), '-struct', 'data');
            stack.geometry.source = "mat";
            testCase.verifyError(@() prepareCOEStack(cfg, stack, temporary.Folder), 'COEReflection:CFCoverage');
        end
    end
end

function cfg = smallConfig()
    cfg = coeReflectionDefaults();
    cfg.periods = [4, 4]; cfg.samplesPerPeriod = 16; cfg.paddingFactor = 2; cfg.outputSize = 129;
    cfg.wavelengths = [450, 550, 650]*1e-9; cfg.spectralWeights = [1, 1, 1]; cfg.channelResponse = eye(3);
end
