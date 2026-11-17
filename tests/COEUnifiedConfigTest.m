classdef COEUnifiedConfigTest < matlab.unittest.TestCase
    methods (TestClassSetup)
        function setup(testCase)
            root=fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'topic_coe_reflection','src')));
            addCoeAlgorithmFixture(testCase,root);
        end
    end
    methods (Test)
        function allActiveProjectsParseAndBuildFourPitch(testCase)
            root=fileparts(fileparts(mfilename('fullpath')));
            inputRoot=fullfile(root,'topic_coe_reflection','Input');
            files=dir(fullfile(inputRoot,'*','config.txt'));
            testCase.verifyEqual(numel(files),25);
            for j=1:numel(files)
                settings=projectSettings();
                [settings,~]=coreLoadProjectConfig(fullfile(files(j).folder,files(j).name),settings);
                cfg=settings.cfg;
                testCase.verifyEqual(cfg.periods,[4,4]); testCase.verifyEqual(cfg.angleHalfRangeDeg,5);
                testCase.verifyEqual(cfg.propagationDistance,0.25); testCase.verifyEqual(cfg.source.type,"plane");
                testCase.verifyFalse(cfg.observer.enabled); testCase.verifyEqual(settings.spectrumMode,"D65");
                testCase.verifyEqual(settings.wavelengthRangeNm,[380,780]); testCase.verifyEqual(settings.wavelengthStepNm,10);
                cfg.wavelengths=(380:10:780)*1e-9;
                [~,geometry]=corePrepareLayerStack(cfg,settings.stack,string(files(j).folder));
                layer=geometry.layerNames(1);
                testCase.verifySize(geometry.layerMasks.(layer),[256,256]);
                testCase.verifyEqual(geometry.selectedPeriods.*geometry.repeatCounts,[4,4]);
            end
        end
        function centralCropPreservesNativePixelGeometry(testCase)
            root=fileparts(fileparts(mfilename('fullpath')));
            folder=fullfile(root,'topic_coe_reflection','Input','Case1');
            settings=projectSettings(); [settings,~]=coreLoadProjectConfig(fullfile(folder,'config.txt'),settings);
            settings.cfg.wavelengths=550e-9;
            [~,geometry]=corePrepareLayerStack(settings.cfg,settings.stack,string(folder));
            image=imread(fullfile(folder,'anode.png'));
            if size(image,3)>1, image=mean(double(image),3); end
            expected=image(385:640,385:640)>0;
            testCase.verifyEqual(geometry.layerMasks.Anode,expected);
            testCase.verifyEqual(geometry.inputPeriods,[16,16]);
            testCase.verifyEqual(geometry.cropBoundsPixels,[385,640,385,640]);
            testCase.verifyEqual(imread(fullfile(folder,'anode.png')),uint8(image));
        end
        function pixelReplicationPreservesOriginalBlocks(testCase)
            root=fileparts(fileparts(mfilename('fullpath')));
            folder=fullfile(root,'topic_coe_reflection','Input','CaseSingleReflectance');
            settings=projectSettings(); [settings,~]=coreLoadProjectConfig(fullfile(folder,'config.txt'),settings);
            settings.cfg.wavelengths=550e-9;
            [~,geometry]=corePrepareLayerStack(settings.cfg,settings.stack,string(folder));
            original=imread(fullfile(folder,'singlelayer.bmp'));
            expected=original>0; if size(expected,3)>1, expected=mean(double(original),3)>0; end
            testCase.verifyEqual(geometry.layerMasks.singlelayer,repelem(expected,4,4));
            testCase.verifyEqual(geometry.pixelReplicationFactor,4);
            testCase.verifyEqual(imread(fullfile(folder,'singlelayer.bmp')),original);
        end
    end
end

function settings=projectSettings()
    stack=struct('layerNames',["Anode","IJP","BM","CF_R","CF_B","CF_G"], ...
        'layerLevels',1:6,'incidentN',1,'exitN',1,'backgroundMode',"absorber");
    stack.geometry=struct('source',"images",'matFile',"layers.mat",'threshold',0.5, ...
        'repeatToFit',false,'inputSamplesPerPeriod',0,'cropToFit',false,'pixelReplication',false);
    for name=stack.layerNames, stack.layers.(name)=coreLayerDefaults(name); end
    settings=struct('cfg',coeReflectionDefaults(),'stack',stack,'spectrumMode',"D65", ...
        'wavelengthRangeNm',[380,780],'wavelengthStepNm',10,'sourceSpectrumFile',"", ...
        'monochromaticWavelengthNm',532,'inputFile',"",'previewWavelengthNm',[620,550,450], ...
        'caseName',"coe_reflection",'saveResults',true,'showProgress',true,'displayFloorDb',-60);
    settings.observation=struct('centerDeg',[0,0],'halfAngleDeg',2.5,'innerHalfAngleDeg',0,'curveSamples',301);
    settings.spotOptions=struct('wavelengthNm',[620,550,450],'zeroOrderExclusionDeg',0.3, ...
        'integrationHalfAngleDeg',0.15,'minSeparationDeg',0.35,'maxSpots',8,'minRelativePeak',0.08);
    settings.displayStyle=struct('showAnnotations',false);
    settings.rgbOptions=struct('mode',"legacy",'clipQuantile',0.9985,'alpha',60,'gamma',0.72,'peakLevel',1);
    settings.colorOptions=struct('toneMapping',"log",'toneStrength',60,'luminanceGamma',0.72,'fixedGain',0,'exposureEV',0);
    settings.colorView=struct('centerDeg',[0,0],'halfAngleDeg',1.5);
end
