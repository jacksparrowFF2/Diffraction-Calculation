classdef COEEyeImagingTest < matlab.unittest.TestCase
    properties (TestParameter)
        sourceQuadrature=struct('cartesian',"cartesian",'polar',"polarGauss");
    end
    methods (TestClassSetup)
        function setup(testCase)
            root=fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'topic_coe_reflection','src')));
            addCoeAlgorithmFixture(testCase,root);
        end
    end
    methods (Test)
        function polarDiskQuadratureMatchesCoherenceAndMoments(testCase)
            [positions,weights]=coreDiskSourceQuadrature(2,7,32);
            testCase.verifyEqual(sum(weights),1,'AbsTol',1e-14);
            testCase.verifyLessThanOrEqual(max(vecnorm(positions,2,2)),1);
            testCase.verifyEqual(weights*positions,[0,0],'AbsTol',1e-14);
            testCase.verifyEqual(weights*(positions(:,1).^2),0.25,'AbsTol',1e-14);
            q=linspace(0,19,401);
            calculated=weights*exp(1i*positions(:,1)*q);
            expected=ones(size(q)); expected(2:end)=2*besselj(1,q(2:end))./q(2:end);
            testCase.verifyLessThan(max(abs(calculated-expected)),1e-4);
        end
        function fourierIntegralPreservesComplexPhase(testCase)
            x=(-3:4)*0.003; y=(-2:3)*0.002;
            fx=linspace(-21,29,11); fy=linspace(-35,17,9);
            [xx,yy]=meshgrid(x,y); field=exp(1i*(17*xx-31*yy)).*(1+4*xx+2*yy);
            expected=exp(-2i*pi*fy(:)*y)*field*exp(-2i*pi*x(:)*fx)*0.003*0.002;
            actual=coreFourierIntegral(field,x,y,fx,fy);
            testCase.verifyEqual(actual,expected,'AbsTol',1e-17);
        end
        function independentSourcesAddIntensities(testCase,sourceQuadrature)
            cfg=eyeConfig(); cfg.source.diameter=0.2e-3; cfg.source.samplesPerAxis=3;
            cfg.source.quadrature=sourceQuadrature; cfg.source.radialNodes=3; cfg.source.azimuthNodes=8;
            mixed=runCOEReflection(cfg); image=mixed.eyeImage;
            expected=zeros(size(image.spectralEnergy));
            positions=image.sourceInfo.sourcePositionsXY;
            single=cfg; single.source.diameter=0;
            for j=1:size(positions,1)
                single.source.centerXY=positions(j,:);
                point=runCOEReflection(single);
                expected=expected+image.sourceWeights(j)*point.eyeImage.spectralEnergy;
            end
            testCase.verifyEqual(image.spectralEnergy,expected,'RelTol',1e-10,'AbsTol',1e-18);
        end
        function phaseAndPowerAreHandledBeforeIntensity(testCase)
            cfg=eyeConfig(); original=runCOEReflection(cfg);
            cfg.complexCoefficient=0.5*exp(0.73i);
            changed=runCOEReflection(cfg);
            testCase.verifyEqual(changed.eyeImage.spectralEnergy,original.eyeImage.spectralEnergy/4, ...
                'RelTol',1e-10,'AbsTol',1e-18);
            testCase.verifyEqual(changed.eyeImage.pupilPowerFraction,original.eyeImage.pupilPowerFraction/4,'RelTol',1e-12);
            testCase.verifyLessThanOrEqual(original.eyeImage.windowToPupilFraction,1.002);
        end
        function focusIsAppliedAtAmplitudeLevel(testCase)
            cfg=eyeConfig(); virtual=runCOEReflection(cfg);
            cfg.observer.focusMode="screen"; screen=runCOEReflection(cfg);
            testCase.verifyEqual(virtual.eyeImage.focusDistance,0.5,'AbsTol',1e-14);
            testCase.verifyEqual(screen.eyeImage.focusDistance,0.25,'AbsTol',1e-14);
            [x,y]=meshgrid(virtual.eyeImage.coordinate); centre=hypot(x,y)<deg2rad(0.03);
            virtualCentre=sum(virtual.eyeImage.spectralEnergy(centre))/sum(virtual.eyeImage.spectralEnergy,'all');
            screenCentre=sum(screen.eyeImage.spectralEnergy(centre))/sum(screen.eyeImage.spectralEnergy,'all');
            testCase.verifyGreaterThan(screenCentre,virtualCentre+0.1);
            testCase.verifyEqual(screen.eyeImage.pupilPowerFraction,virtual.eyeImage.pupilPowerFraction,'RelTol',1e-13);
        end
        function rejectedSamplingDoesNotSilentlyAlias(testCase)
            cfg=eyeConfig(); cfg.observer.angleHalfRangeDeg=4.5;
            testCase.verifyError(@()runCOEReflection(cfg),'COEReflection:EyePupilSampling');
            cfg.observer.pupilSamples=1537;
            testCase.verifyError(@()runCOEReflection(cfg),'COEReflection:EyeImageSampling');
        end
        function virtualMirrorSourceHasCorrectImageSignAndScale(testCase)
            cfg=eyeConfig(); cfg.periods=[80,80]; cfg.samplesPerPeriod=8;
            cfg.source.centerXY=[0.3e-3,0];
            mirror=runCOEReflection(cfg); view=mirror.eyeImage;
            [~,index]=max(view.spectralEnergy,[],'all');
            [row,column]=ind2sub(size(view.spectralEnergy),index);
            testCase.verifyEqual(view.coordinate(column),-0.3e-3/0.5,'AbsTol',deg2rad(0.003));
            testCase.verifyEqual(view.coordinate(row),0,'AbsTol',deg2rad(0.0015));
        end
        function constantHeightChangesPhaseNotIntensity(testCase)
            cfg=eyeConfig(); cfg.heightMap=25e-9; result=runCOEReflection(cfg);
            testCase.verifyTrue(result.eyeImage.heightPhaseEnabled);
            testCase.verifyEqual(result.eyeImage.spectralEnergy,result.eyeImage.baselineSpectralEnergy, ...
                'RelTol',1e-10,'AbsTol',1e-18);
        end
        function fixedDisplayReferencePreservesIntensityScaling(testCase)
            signal=reshape(linspace(0,2,75),5,5,3);
            options=struct('mode',"legacy",'clipQuantile',0.9,'alpha',60,'gamma',0.72,'fixedReferenceSignal',4);
            [small,~,mapping]=plotComposeReflectionRgb(signal,signal,options,false);
            padded=zeros(15,15,3); padded(6:10,6:10,:)=signal;
            [large,~,largeMapping]=plotComposeReflectionRgb(padded,padded,options,false);
            testCase.verifyEqual(small,large(6:10,6:10,:),'AbsTol',1e-14);
            testCase.verifyEqual(mapping.referenceSignal,largeMapping.referenceSignal);
            testCase.verifyEqual(mapping.referenceSignal,4);
        end
        function legacyConfigurationWithoutObserverStillRuns(testCase)
            cfg=rmfield(eyeConfig(),'observer'); cfg.angleHalfRangeDeg=0.15; cfg.metricAngleDeg=0.1;
            result=runCOEReflection(cfg);
            testCase.verifyTrue(isfield(result,'spectralEnergy'));
            testCase.verifyFalse(isfield(result,'eyeImage'));
        end
        function blackReflectorHasNoSpuriousImage(testCase)
            cfg=eyeConfig(); cfg.complexCoefficient=0; result=runCOEReflection(cfg);
            testCase.verifyEqual(result.eyeImage.spectralEnergy,zeros(size(result.eyeImage.spectralEnergy)));
            testCase.verifyEqual(result.eyeImage.pupilPowerFraction,0);
            testCase.verifyTrue(all(isfinite(result.eyeImage.windowToPupilFraction)));
        end
    end
end

function cfg=eyeConfig()
    cfg=coeReflectionDefaults(); cfg.periods=[4,4]; cfg.samplesPerPeriod=16;
    cfg.complexCoefficient=1; cfg.wavelengths=550e-9; cfg.spectralWeights=1; cfg.channelResponse=1;
    cfg.source.type="point"; cfg.source.wavefront="paraxial"; cfg.source.distance=0.25;
    cfg.source.diameter=0; cfg.observer.enabled=true; cfg.observer.model="coherentEye";
    cfg.observer.angleHalfRangeDeg=0.15; cfg.outputSize=201;
    cfg.observer.pupilSamples=129;
end
