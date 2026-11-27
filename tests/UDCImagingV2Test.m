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
        function specifiedTransmittanceAppliesToWhitePixels(testCase)
            stack = testCase.stackFor("specifiedTransmittance");
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
            index = 1.5+0.2i;
            r = (1-index)/(1+index);
            phase = exp(2i*pi*index*film.thickness./testCase.Config.wavelengths);
            t = (1-r^2)*phase./(1-r^2*phase.^2);
            testCase.verifyEqual(optics.transmissionLookup(optics.regionPatterns,:),t,AbsTol=1e-13);
            testCase.verifyEqual(optics.transmissionLookup(~optics.regionPatterns,:),zeros(1,3),AbsTol=1e-14);
            testCase.verifyEqual(optics.mtf(17,17,:),ones(1,1,3),AbsTol=1e-13);
        end

        function nkFileInterpolatesAndMatchesInlineIndices(testCase)
            stack = testCase.stackFor("nk");
            temporary = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture());
            folder = string(temporary.Folder);
            copyfile(fullfile(testCase.InputDir,stack.layers.Opening.file),folder);
            writematrix([400,1.4,0.1;700,1.7,0.4],fullfile(folder,'spectrum.csv'));
            stack.layers.Opening.nkFile = "spectrum.csv";
            fromFile = runUDCImagingV2(testCase.Config,stack,folder);
            stack.layers.Opening.nkFile = "";
            stack.layers.Opening.n = [1.45,1.55,1.65];
            stack.layers.Opening.k = [0.15,0.25,0.35];
            fromInline = runUDCImagingV2(testCase.Config,stack,folder);
            testCase.verifyEqual(fromFile.transmissionLookup,fromInline.transmissionLookup,AbsTol=1e-13);
        end

        function nkFileRejectsWavelengthOutsideTable(testCase)
            cfg = testCase.Config;
            cfg.wavelengths(end) = 790e-9;
            stack = testCase.stackFor("nk");
            testCase.verifyError(@() runUDCImagingV2(cfg,stack,testCase.InputDir),'UDCImagingV2:NKRange');
        end

        function croppedFftMatchesFullFftBinning(testCase)
            cfg = testCase.Config;
            cfg.wavelengths = 550e-9;
            cfg.spectralWeights = 1;
            cfg.channelResponse = 1;
            cfg.receiverMethod = "fftBins";
            cfg.receiverFftPaddingFactor = 2;
            optics = runUDCImagingV2(cfg,testCase.stackFor("nk"),testCase.InputDir);
            coefficient = reshape(optics.transmissionLookup(double(optics.regionMap(:))+1),size(optics.pupil));
            field = optics.pupil.*coefficient;
            nfft = max(size(field))*cfg.receiverFftPaddingFactor;
            spectrum = fftshift(fft2(single(field),nfft,nfft));
            native = (-nfft/2:nfft/2-1)*cfg.wavelengths*cfg.focalLength/(nfft*optics.layerGeometry.gridStep);
            step = native(2)-native(1);
            distance = abs(optics.coordinate(:)-native);
            width = min(cfg.sensorPitch,step);
            weights = max(0,min(width,(cfg.sensorPitch+step)/2-distance))/step;
            reference = weights*double(abs(spectrum).^2)*weights.'/(nfft^2*nnz(optics.pupil));
            testCase.verifyEqual(optics.spectralEnergy,reference,AbsTol=1e-8);
        end

        function focusedObjectProjectionUsesDistanceAndSensorPitch(testCase)
            cfg = testCase.Config;
            cfg.receiverMethod = "fftBins";
            cfg.receiverFftPaddingFactor = 2;
            cfg.sensorPitch = 1e-6;
            cfg.sceneGeometry = struct('mode',"object",'objectSizeMm',[0.9,0.45],'objectDistanceMm',30);
            scene = ones(20,40);
            stack = testCase.stackFor("nk");
            [near,nearImages] = runUDCImagingV2(cfg,stack,testCase.InputDir,scene);
            testCase.verifySize(nearImages.inputLinear,[50,100,3]);
            testCase.verifyEqual(near.config.sceneProjection.imageDistanceMm,10/3,AbsTol=1e-12);
            testCase.verifyEqual(near.config.sceneProjection.magnification,1/9,AbsTol=1e-12);
            cfg.sceneGeometry.objectDistanceMm = 57;
            [~,farImages] = runUDCImagingV2(cfg,stack,testCase.InputDir,scene);
            testCase.verifySize(farImages.inputLinear,[25,50,3]);
            cfg.sceneGeometry.mode = "sensor";
            [sensor,sensorImages] = runUDCImagingV2(cfg,stack,testCase.InputDir,scene);
            testCase.verifySize(sensorImages.inputLinear,[20,40,3]);
            testCase.verifyNotEqual(sensor.psf,near.psf);
        end

        function sharedFftPreservesEachScenesBinning(testCase)
            cfg = testCase.Config;
            cfg.receiverMethod = "fftBins";
            cfg.receiverFftPaddingFactor = 2;
            cfg.sensorPitch = 1e-6;
            cfg.sceneGeometry = struct('mode',"object",'objectSizeMm',[0.9,0.45],'objectDistanceMm',30);
            scene = reshape(linspace(0,1,800),20,40);
            stack = testCase.stackFor("nk");
            scenes = ["object","sensor"];
            [shared,sharedImages] = runUDCImagingV2(cfg,stack,testCase.InputDir,scene,"linear",scenes);
            testCase.verifySize(shared,[1,2]);
            for k = 1:numel(scenes)
                cfg.sceneGeometry.mode = scenes(k);
                [single,singleImages] = runUDCImagingV2(cfg,stack,testCase.InputDir,scene);
                for field = ["spectralEnergy","clearKernel","psf","mtf"]
                    testCase.verifyEqual(shared(k).(field),single.(field),AbsTol=1e-10);
                end
                testCase.verifyEqual(sharedImages(k).fixedExposure,singleImages.fixedExposure,AbsTol=1e-10);
                testCase.verifyEqual(sharedImages(k).brightnessMatched,singleImages.brightnessMatched,AbsTol=1e-10);
            end
        end

        function framedPreviewMatchesNativeConvolutionAndAveraging(testCase)
            cfg = testCase.Config;
            cfg.receiverMethod = "fftBins";
            cfg.receiverFftPaddingFactor = 2;
            cfg.sensorPitch = 1e-6;
            cfg.sceneGeometry = struct('mode',"object",'objectSizeMm',[0.9,0.45],'objectDistanceMm',30, ...
                'sensorFrameMm',[0.16,0.10],'previewWidthPixels',80);
            [optics,images] = runUDCImagingV2(cfg,testCase.stackFor("nk"),testCase.InputDir,ones(20,40));
            testCase.verifySize(images.inputLinear,[50,80,3]);
            testCase.verifyEqual(optics.config.sceneProjection.frame.sizePixels,[160,100]);
            background = 0.18;
            target = zeros(100,160);
            target(26:75,31:130) = 1-background;
            native = background+conv2(target,optics.psf(:,:,2),'same');
            expected = (native(1:2:end,1:2:end)+native(2:2:end,1:2:end)+ ...
                native(1:2:end,2:2:end)+native(2:2:end,2:2:end))/4;
            testCase.verifyEqual(images.brightnessMatched(:,:,2),expected,AbsTol=1e-12);
            testCase.verifyEqual(images.inputLinear(1,1,2),background,AbsTol=1e-12);
        end

        function rejectsGrayBeforeResampling(testCase)
            stack = testCase.stackFor("specifiedTransmittance");
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
