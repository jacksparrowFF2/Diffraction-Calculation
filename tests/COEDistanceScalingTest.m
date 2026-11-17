classdef COEDistanceScalingTest < matlab.unittest.TestCase
    %COEDistanceScalingTest Public-API peak positions, energy mapping and limits.
    properties (TestParameter)
        distances=struct('nearSource',[0.15,0.25,0.375], ...
            'equalDistances',[0.25,0.25,0.5], ...
            'farSource',[0.35,0.25,7/12], ...
            'nearObserver',[0.25,0.15,0.625], ...
            'farObserver',[0.25,0.35,5/12]);
    end
    methods (TestClassSetup)
        function addModule(testCase)
            root=fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root,'topic_coe_reflection','src')));
        end
    end
    methods (Test)
        function phaseRampPeakFollowsDistances(testCase,distances)
            cfg=phaseRampConfig();
            cfg.source.distance=distances(1);
            cfg.observer.distance=distances(2);
            result=runCOEReflection(cfg);
            view=result.distanceView;
            [~,linear]=max(view.spectralEnergy(:));
            [row,col]=ind2sub(size(view.spectralEnergy),linear);
            step=view.coordinate(2)-view.coordinate(1);

            testCase.verifyEqual(view.angularScale,distances(3),'AbsTol',1e-14);
            testCase.verifyEqual(view.coordinate(row),0,'AbsTol',step/2);
            testCase.verifyEqual(view.coordinate(col), ...
                distances(3)*550e-9/50e-6,'AbsTol',step/2);
            testCase.verifyEqual(view.screenCoordinate(col), ...
                distances(2)*distances(3)*550e-9/50e-6,'AbsTol',distances(2)*step/2);
        end
        function energyBinsAndDensityUseMappedArea(testCase)
            cfg=phaseRampConfig();
            result=runCOEReflection(cfg);
            view=result.distanceView;
            reference=cfg;
            reference.observer.enabled=false;
            reference.source.type="plane";
            reference.angleHalfRangeDeg=2;
            plane=runCOEReflection(reference);
            referenceBin=(plane.coordinate(2)-plane.coordinate(1))^2;

            testCase.verifyEqual(view.spectralEnergy,plane.spectralEnergy,'AbsTol',1e-13);
            testCase.verifyEqual(view.spectralAngularDensity, ...
                plane.spectralEnergy/referenceBin/0.5^2,'RelTol',1e-12,'AbsTol',1e-9);
            testCase.verifyEqual(view.referenceThroughput,1,'AbsTol',1e-13);
            testCase.verifyEqual(result.metrics.Throughput,plane.metrics.Throughput,'AbsTol',1e-13);
            testCase.verifyEqual(view.referencePropagationRuns,1);
        end
        function narrowWindowDoesNotRenormalizeLostPeak(testCase)
            cfg=phaseRampConfig();
            cfg.observer.angleHalfRangeDeg=0.1;
            result=runCOEReflection(cfg);

            testCase.verifyLessThan(sum(result.distanceView.spectralEnergy,'all'),0.01);
            testCase.verifyEqual(result.distanceView.referenceThroughput,1,'AbsTol',1e-13);
        end
        function distanceViewDoesNotChangeFarField(testCase)
            cfg=phaseRampConfig();
            enabled=runCOEReflection(cfg);
            cfg.observer.enabled=false;
            disabled=runCOEReflection(cfg);

            testCase.verifyEqual(enabled.spectralEnergy,disabled.spectralEnergy,'AbsTol',1e-13);
            testCase.verifyEqual(enabled.spectralIrradiance,disabled.spectralIrradiance,'AbsTol',1e-13);
            testCase.verifyFalse(isfield(disabled,'distanceView'));
        end
        function offAxisSourceIsRejected(testCase)
            cfg=phaseRampConfig(); cfg.source.centerXY=[0.001,0];
            testCase.verifyError(@() runCOEReflection(cfg),'COEReflection:DistanceViewSource');
        end
        function finiteSourceMatchesDirectLinearConvolution(testCase)
            cfg=phaseRampConfig(); cfg.source.diameter=1e-3;
            result=runCOEReflection(cfg); view=result.distanceView;
            plane=runCOEReflection(view.referenceConfig);
            expected=conv2(plane.spectralEnergy,view.kernels.sourceKernel,'same');
            selected=view.kernels.paddingBins+(1:cfg.outputSize);

            testCase.verifyEqual(view.spectralEnergy,expected(selected,selected),'AbsTol',1e-13);
            testCase.verifyEqual(sum(view.kernels.sourceKernel,'all'),1,'AbsTol',1e-14);
            testCase.verifyEqual(view.kernels.sourceAngularDiameterRad,0.002,'AbsTol',1e-14);
            testCase.verifyLessThan(max(view.convolutionMassError,[],'all'),1e-12);
        end
        function circularSourceHasExpectedSecondMoment(testCase)
            cfg=phaseRampConfig(); cfg.source.diameter=1e-3;
            result=runCOEReflection(cfg); kernels=result.distanceView.kernels;
            half=(size(kernels.sourceKernel,1)-1)/2;
            [x,y]=meshgrid((-half:half)*kernels.stepRad);
            secondMoment=sum((x.^2+y.^2).*kernels.sourceKernel,'all');

            testCase.verifyEqual(secondMoment,(0.001)^2/2,'AbsTol',2*kernels.stepRad^2);
        end
        function sourceAveragingIncludesPointPeakOutsideView(testCase)
            cfg=phaseRampConfig(); cfg.source.diameter=5e-3;
            cfg.observer.angleHalfRangeDeg=0.2;
            result=runCOEReflection(cfg); view=result.distanceView;

            testCase.verifyLessThan(sum(view.pointSpectralEnergy,'all'),0.01);
            testCase.verifyGreaterThan(sum(view.spectralEnergy,'all'),0.05);
            testCase.verifyGreaterThan(view.kernels.paddingBins,0);
        end
        function mirrorRegionUsesSourceAndPupilIndependently(testCase,distances)
            cfg=phaseRampConfig(); cfg.source.diameter=5e-3;
            cfg.source.distance=distances(1); cfg.observer.distance=distances(2);
            result=runCOEReflection(cfg); region=result.distanceView.participation;
            expected=(distances(2)*5e-3+distances(1)*3e-3)/(distances(1)+distances(2));

            testCase.verifyEqual(region.diameter,expected,'AbsTol',1e-14);
            testCase.verifyEqual(region.nominalPeriods,expected/cfg.pitch,'AbsTol',1e-12);
            testCase.verifyEqual(region.geometricArea,pi*(expected/2)^2,'AbsTol',1e-16);
            testCase.verifyEqual(sum(region.rayPairDensity,'all')*region.diagnosticStep^2,1,'AbsTol',3e-4);
            testCase.verifyLessThan(region.rayWeightInModel,1);
            testCase.verifyLessThan(region.modeledArea,region.geometricArea);
            testCase.verifyEqual(result.metrics.Throughput,1,'AbsTol',1e-12);
        end
        function airyKernelRetainsAnalyticMassAndDoesNotNormalizeTail(testCase)
            cfg=phaseRampConfig(); cfg.observer.pupilModel="airyIntensity";
            result=runCOEReflection(cfg); kernels=result.distanceView.kernels;
            z=cfg.observer.airyRadiusFactor*3.8317059702075125;
            mass=1-besselj(0,z)^2-besselj(1,z)^2;

            testCase.verifyEqual(sum(kernels.pupilKernels{1},'all'),mass,'AbsTol',1e-13);
            testCase.verifyLessThan(kernels.pupilRetainedMass,1);
            testCase.verifyEqual(kernels.airyFirstZeroRad, ...
                1.2196698912665045*550e-9/3e-3,'AbsTol',1e-14);
            testCase.verifyGreaterThan(result.distanceView.kernels.paddingBins,0);
        end
        function largerPupilNarrowsAiryAndExpandsGeometricRegion(testCase)
            cfg=phaseRampConfig(); cfg.observer.pupilModel="airyIntensity";
            small=runCOEReflection(cfg);
            cfg.observer.pupilDiameter=6e-3;
            large=runCOEReflection(cfg);

            testCase.verifyEqual(large.distanceView.kernels.airyFirstZeroRad, ...
                small.distanceView.kernels.airyFirstZeroRad/2,'AbsTol',1e-14);
            testCase.verifyEqual(large.distanceView.participation.diameter, ...
                2*small.distanceView.participation.diameter,'AbsTol',1e-14);
            testCase.verifyEqual(large.spectralEnergy,small.spectralEnergy,'AbsTol',1e-13);
        end
        function zeroPupilAndSourceRecoverLegacyPointView(testCase)
            cfg=phaseRampConfig(); cfg.observer.pupilDiameter=0;
            cfg.observer.pupilModel="airyIntensity";
            result=runCOEReflection(cfg); plane=runCOEReflection(result.distanceView.referenceConfig);

            testCase.verifyEqual(result.distanceView.spectralEnergy,plane.spectralEnergy,'AbsTol',1e-13);
            testCase.verifyEqual(result.distanceView.participation.diameter,0,'AbsTol',1e-14);
        end
        function fftBinsApproximationPreservesBroadSourcePattern(testCase)
            cfg=phaseRampConfig(); cfg.source.diameter=5e-3;
            quadrature=runCOEReflection(cfg);
            cfg.observer.referenceMethod="fftBins";
            cfg.observer.referencePaddingFactor=2;
            fftResult=runCOEReflection(cfg);
            exact=quadrature.distanceView.spectralEnergy;
            approximate=fftResult.distanceView.spectralEnergy;
            relativeL1=sum(abs(approximate-exact),'all')/sum(exact,'all');

            testCase.verifyLessThan(relativeL1,0.04);
            testCase.verifyEqual(fftResult.distanceView.referenceThroughput,1,'AbsTol',1e-12);
            testCase.verifyEqual(fftResult.spectralEnergy,quadrature.spectralEnergy,'AbsTol',1e-13);
        end
        function excessiveReferenceAngleIsRejected(testCase)
            cfg=phaseRampConfig(); cfg.observer.distance=3;
            testCase.verifyError(@() runCOEReflection(cfg),'COEReflection:DistanceViewRange');
        end
        function negativeObserverDistanceIsRejected(testCase)
            cfg=phaseRampConfig(); cfg.observer.distance=-0.25;
            testCase.verifyError(@() runCOEReflection(cfg),'MATLAB:expectedPositive');
        end
        function copiedModuleResolvesDistanceView(testCase)
            temporary=testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture());
            root=fileparts(fileparts(mfilename('fullpath')));
            folder=fullfile(temporary.Folder,'coe'); mkdir(folder);
            copyfile(fullfile(root,'topic_coe_reflection','src'),fullfile(folder,'src'));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(folder,'src')));
            cfg=phaseRampConfig();
            result=runCOEReflection(cfg);

            testCase.verifyEqual(string(which('runCOEReflection')), ...
                string(fullfile(folder,'src','runCOEReflection.m')));
            testCase.verifyEqual(result.distanceView.angularScale,0.5,'AbsTol',1e-14);
        end
    end
end

function cfg=phaseRampConfig()
    cfg=coeReflectionDefaults();
    cfg.periods=[16,16]; cfg.samplesPerPeriod=16;
    cfg.wavelengths=550e-9; cfg.spectralWeights=1; cfg.channelResponse=1;
    cfg.outputSize=401; cfg.angleHalfRangeDeg=2; cfg.metricAngleDeg=1;
    count=cfg.periods(1)*cfg.samplesPerPeriod;
    x=(0:count-1)*cfg.pitch/cfg.samplesPerPeriod;
    cfg.complexCoefficient=repmat(exp(2i*pi*x/cfg.pitch),count,1);
    cfg.source.type="point"; cfg.source.distance=0.25;
    cfg.source.diameter=0; cfg.source.integrationMethod="direct";
    cfg.source.apertureShape="none";
    cfg.observer.enabled=true; cfg.observer.distance=0.25;
    cfg.observer.angleHalfRangeDeg=1;
end
