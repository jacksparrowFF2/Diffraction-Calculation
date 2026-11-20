classdef UDCImagingV2Test < matlab.unittest.TestCase
    properties (Access=private)
        InputDir
        Settings
        Config
    end

    methods (TestClassSetup)
        function loadCase01(testCase)
            root = fileparts(fileparts(mfilename('fullpath')));
            module = fullfile(root,'topic_udc_imagingV2');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(module,'src')));
            testCase.InputDir = string(fullfile(module,'Input','case01'));
            testCase.Settings = jsondecode(fileread(fullfile(testCase.InputDir,'config.json')));
            cfg = testCase.Settings.cfg;
            cfg.periods = reshape(cfg.periods,1,[]);
            cfg.wavelengths = [450,550,650]*1e-9;
            cfg.spectralWeights = [1,1,1];
            cfg.channelResponse = eye(3);
            cfg.outputSize = 33;
            testCase.Config = cfg;
        end
    end

    methods (Test)
        function directAppliesPowerTransmittanceToWhitePixels(testCase)
            stack = testCase.stackFor("direct");
            optics = runUDCImagingV2(testCase.Config,stack,testCase.InputDir);
            openFraction = nnz(optics.layerGeometry.unionMask & optics.pupil)/nnz(optics.pupil);
            testCase.verifyGreaterThan(openFraction,0);
            testCase.verifyLessThan(openFraction,1);
            testCase.verifyEqual(optics.transmissionLookup(~optics.regionPatterns,:),zeros(1,3),AbsTol=1e-14);
            testCase.verifyEqual(optics.transmissionLookup(optics.regionPatterns,:), ...
                sqrt(stack.layers.Opening.transmittance)*ones(1,3),AbsTol=1e-14);
            testCase.verifyEqual(optics.metrics.Throughput, ...
                openFraction*stack.layers.Opening.transmittance*ones(3,1),AbsTol=1e-13);
        end

        function nkMatchesComplexSingleFilmSolution(testCase)
            stack = testCase.stackFor("nk");
            optics = runUDCImagingV2(testCase.Config,stack,testCase.InputDir);
            film = stack.layers.Opening;
            index = film.n+1i*film.k;
            r = (1-index)/(1+index);
            phase = exp(2i*pi*index*film.thickness./testCase.Config.wavelengths);
            t = (1-r^2)*phase./(1-r^2*phase.^2);
            testCase.verifyEqual(optics.transmissionLookup(optics.regionPatterns,:),t,AbsTol=1e-13);
            testCase.verifyEqual(optics.transmissionLookup(~optics.regionPatterns,:),zeros(1,3),AbsTol=1e-14);
            testCase.verifyEqual(optics.mtf(17,17,:),ones(1,1,3),AbsTol=1e-13);
        end

        function rejectsGrayBeforeResampling(testCase)
            stack = testCase.stackFor("direct");
            temporary = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture());
            mask = imread(fullfile(testCase.InputDir,stack.layers.Opening.file));
            mask(1,1) = 128; % 原图第一像素会在配置的降采样中被跳过。
            imwrite(mask,fullfile(temporary.Folder,stack.layers.Opening.file));
            testCase.verifyError(@() runUDCImagingV2(testCase.Config,stack,string(temporary.Folder)), ...
                'UDCImagingV2:BinaryMask');
        end
    end

    methods (Access=private)
        function stack = stackFor(testCase, mode)
            stack = testCase.Settings.stack;
            stack.layerNames = reshape(string(stack.layerNames),1,[]);
            stack.layerLevels = reshape(stack.layerLevels,1,[]);
            stack.layers = testCase.Settings.transmissionModes.(mode);
        end
    end
end
