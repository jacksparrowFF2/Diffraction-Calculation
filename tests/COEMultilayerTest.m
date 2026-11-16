classdef COEMultilayerTest < matlab.unittest.TestCase
    %COEMultilayerTest Analytic film limits, absorber footprint and layer input contracts.
    methods (TestClassSetup)
        function addModule(testCase)
            root = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'topic_coe_reflection', 'src')));
            addCoeAlgorithmFixture(testCase, root);
        end
    end
    methods (Test)
        function rectanglePlaneIlluminationUsesLitAreaAndKnownReflectance(testCase)
            cfg=smallConfig(); cfg.wavelengths=550e-9; cfg.spectralWeights=1; cfg.channelResponse=1;
            cfg.complexCoefficient=0.5; cfg.incidentIrradiance=7;
            cfg.source.apertureShape="rectangle";
            cfg.source.rectangleSize=[100e-6,75e-6];
            value=runCOEReflection(cfg);
            [xx,yy]=meshgrid(value.x,value.y);
            expected=xx>=-50e-6 & xx<50e-6 & yy>=-37.5e-6 & yy<37.5e-6;
            testCase.verifyEqual(value.pupil,double(expected));
            testCase.verifyEqual(value.incidentIrradianceMap,7*double(expected),'AbsTol',1e-12);
            testCase.verifyEqual(value.incidentArea,nnz(expected)*(cfg.pitch/cfg.samplesPerPeriod)^2,'AbsTol',1e-18);
            testCase.verifyEqual(value.incidentPower,7*value.incidentArea,'AbsTol',1e-18);
            testCase.verifyEqual(value.metrics.Throughput,0.25,'AbsTol',1e-12);
            testCase.verifyLessThan(value.metrics.ParsevalAbsoluteError,1e-12);
        end

        function circlePlaneIlluminationUsesDiameterAndCentre(testCase)
            cfg=smallConfig(); cfg.complexCoefficient=0.5;
            cfg.source.apertureShape="circle"; cfg.source.circleDiameter=100e-6;
            cfg.source.centerXY=[12.5e-6,-6.25e-6];
            value=runCOEReflection(cfg);
            [xx,yy]=meshgrid(value.x,value.y);
            expected=(xx-12.5e-6).^2+(yy+6.25e-6).^2<=(50e-6)^2;
            testCase.verifyEqual(value.pupil,double(expected));
            testCase.verifyEqual(value.metrics.Throughput,0.25*ones(3,1),'AbsTol',1e-12);
            testCase.verifyLessThan(abs(value.incidentArea-pi*(50e-6)^2), ...
                2*pi*50e-6*cfg.pitch/cfg.samplesPerPeriod);
        end

        function millimetreBeamDoesNotExpandSmallInput(testCase)
            cfg=smallConfig(); cfg.complexCoefficient=0.5;
            reference=runCOEReflection(cfg);
            cfg.source.apertureShape="rectangle"; cfg.source.rectangleSize=[1e-3,2e-3];
            value=runCOEReflection(cfg);
            testCase.verifyTrue(value.illuminationAperture.extendsBeyondGrid);
            testCase.verifyEqual(value.illuminationAperture.requestedArea,2e-6,'AbsTol',1e-18);
            testCase.verifyEqual(value.incidentArea,reference.incidentArea,'AbsTol',1e-18);
            testCase.verifyEqual(value.spectralEnergy,reference.spectralEnergy,'AbsTol',1e-12);
            testCase.verifyEqual(value.config.periods,cfg.periods);
        end

        function emptyOrInvalidPlaneFootprintIsRejected(testCase)
            cfg=smallConfig(); cfg.source.apertureShape="circle";
            cfg.source.circleDiameter=10e-6; cfg.source.centerXY=[1e-3,0];
            testCase.verifyError(@()runCOEReflection(cfg),'COEReflection:EmptyIllumination');
            cfg.source.apertureShape="triangle";
            testCase.verifyError(@()runCOEReflection(cfg),'MATLAB:validators:mustBeMember');
        end

        function sourceConvolutionMatchesDirectParaxialDisk(testCase)
            [cfg,stack,folder]=sourceConvolutionInput(testCase);
            direct=runCOEReflection(cfg,stack,folder);
            cfg.source.integrationMethod="angularConvolution";
            convolved=runCOEReflection(cfg,stack,folder);
            testCase.verifyEqual(convolved.spectralEnergy,direct.spectralEnergy,'AbsTol',1e-11);
            testCase.verifyEqual(convolved.clearKernel,direct.clearKernel,'AbsTol',1e-11);
            testCase.verifyEqual(convolved.metrics.Throughput,direct.metrics.Throughput,'AbsTol',1e-12);
            testCase.verifyEqual(convolved.sourceConvolution.kernelMass,1,'AbsTol',1e-14);
            testCase.verifyLessThan(convolved.sourceConvolution.massError,1e-12);
        end

        function sourceConvolutionMatchesTiltedPlaneDisk(testCase)
            [cfg,stack,folder]=sourceConvolutionInput(testCase);
            cfg.source.wavefront="plane";
            direct=runCOEReflection(cfg,stack,folder);
            cfg.source.integrationMethod="angularConvolution";
            convolved=runCOEReflection(cfg,stack,folder);
            testCase.verifyEqual(convolved.spectralEnergy,direct.spectralEnergy,'AbsTol',1e-11);
        end

        function sourceConvolutionPadsOutsideObservationWindow(testCase)
            [cfg,stack,folder]=sourceConvolutionInput(testCase);
            cfg.source.centerXY=[deg2rad(1.4)*cfg.source.distance,0];
            direct=runCOEReflection(cfg,stack,folder);
            cfg.source.integrationMethod="angularConvolution";
            convolved=runCOEReflection(cfg,stack,folder);
            testCase.verifyEqual(convolved.spectralEnergy,direct.spectralEnergy,'AbsTol',1e-11);
            testCase.verifyGreaterThan(convolved.sourceConvolution.paddedHalfAngleDeg,cfg.angleHalfRangeDeg);
            testCase.verifyLessThan(convolved.metrics.CapturedThroughput,convolved.metrics.Throughput);
        end

        function zeroDiameterConvolutionReturnsPointSource(testCase)
            [cfg,stack,folder]=sourceConvolutionInput(testCase);
            cfg.source.diameter=0;
            direct=runCOEReflection(cfg,stack,folder);
            cfg.source.integrationMethod="angularConvolution";
            convolved=runCOEReflection(cfg,stack,folder);
            testCase.verifyEqual(convolved.spectralEnergy,direct.spectralEnergy,'AbsTol',1e-11);
            testCase.verifyEqual(convolved.sourceConvolution.kernel,1);
        end

        function exactSphericalConvolutionRequiresExplicitApproximation(testCase)
            [cfg,stack,folder]=sourceConvolutionInput(testCase);
            cfg.source.integrationMethod="angularConvolution";
            cfg.source.wavefront="spherical";
            testCase.verifyError(@()runCOEReflection(cfg,stack,folder),'COEReflection:SourceConvolution');
        end

        function explicitRasterSamplingKeepsPhysicalExtentAndInputFiles(testCase)
            [cfg,stack,folder]=layerInput(testCase);
            original=imread(fullfile(folder,'CF_R.png'));
            stack.geometry.inputSamplesPerPeriod=16;
            cfg.samplesPerPeriod=8; cfg.angleHalfRangeDeg=1.5; cfg.metricAngleDeg=1;
            value=runCOEReflection(cfg,stack,folder);
            expected=original(2:2:end,2:2:end)>0;
            testCase.verifyEqual(value.layerGeometry.layerMasks.CF_R,expected);
            testCase.verifyEqual(value.layerGeometry.inputGridSize,[64,64]);
            testCase.verifyEqual(value.layerGeometry.sampledInputGridSize,[32,32]);
            testCase.verifyEqual(value.layerGeometry.inputPeriods,[4,4]);
            testCase.verifyEqual(value.layerGeometry.gridStep,cfg.pitch/8,'AbsTol',1e-18);
            testCase.verifyEqual(imread(fullfile(folder,'CF_R.png')),original);
        end

        function explicitRasterSamplingCanRepeatCompleteCell(testCase)
            [cfg,stack,folder]=layerInput(testCase);
            stack.geometry.inputSamplesPerPeriod=16;
            stack.geometry.repeatToFit=true;
            cfg.samplesPerPeriod=8; cfg.angleHalfRangeDeg=1.5; cfg.metricAngleDeg=1; cfg.periods=[8,8];
            value=runCOEReflection(cfg,stack,folder);
            testCase.verifyEqual(value.layerGeometry.repeatCounts,[2,2]);
            testCase.verifySize(value.layerGeometry.layerMasks.IJP,[64,64]);
            testCase.verifyTrue(all(value.layerGeometry.layerMasks.IJP,'all'));
        end

        function explicitRasterSamplingRejectsUpsampling(testCase)
            [cfg,stack,folder]=layerInput(testCase);
            stack.geometry.inputSamplesPerPeriod=16;
            cfg.samplesPerPeriod=32;
            testCase.verifyError(@()runCOEReflection(cfg,stack,folder),'COEReflection:RasterUpsampling');
        end

        function repeatTwoPitchMatchesCompleteImage(testCase)
            [repeated,complete]=repeatInputPair(testCase,[2,2]);
            testCase.verifyEqual(repeated.spectralEnergy,complete.spectralEnergy,'AbsTol',1e-12);
            testCase.verifyEqual(repeated.layerGeometry.repeatCounts,[4,4]);
            testCase.verifyEqual(repeated.layerGeometry.inputPeriods,[2,2],'AbsTol',1e-12);
        end

        function repeatFourPitchMatchesCompleteImage(testCase)
            [repeated,complete]=repeatInputPair(testCase,[4,4]);
            testCase.verifyEqual(repeated.spectralEnergy,complete.spectralEnergy,'AbsTol',1e-12);
            testCase.verifyEqual(repeated.layerGeometry.repeatCounts,[2,2]);
        end

        function rectangularUnitUsesSeparateRepeatCounts(testCase)
            [repeated,complete]=repeatInputPair(testCase,[2,4]);
            testCase.verifyEqual(repeated.spectralEnergy,complete.spectralEnergy,'AbsTol',1e-12);
            testCase.verifyEqual(repeated.layerGeometry.repeatCounts,[4,2]);
        end

        function nonDivisibleTargetIsRejected(testCase)
            [cfg,stack,folder]=singleLayerInput(testCase);
            cfg.periods=[6,6]; stack.geometry.repeatToFit=true;
            testCase.verifyError(@()runCOEReflection(cfg,stack,folder),'COEReflection:RepeatExtent');
        end

        function differentLayerUnitSizesAreRejected(testCase)
            [cfg,stack,folder]=layerInput(testCase);
            stack.geometry.repeatToFit=true;
            imwrite(true(32),fullfile(folder,'IJP.png'));
            testCase.verifyError(@()runCOEReflection(cfg,stack,folder),'COEReflection:InputAlignment');
        end

        function extendedSourceEqualsIncoherentWeightedPointSum(testCase)
            [cfg,stack,folder]=singleLayerInput(testCase);
            cfg.wavelengths=550e-9;cfg.spectralWeights=1;cfg.channelResponse=1;
            stack.layers.singlelayer.method="reflectance";
            stack.layers.singlelayer.reflectance=0.36;
            cfg.source.type="point";cfg.source.distance=0.01;cfg.source.diameter=0.2e-3;
            result=runCOEReflection(cfg,stack,folder);
            expected=independentSourcePointSum(cfg,stack,folder);
            testCase.verifyEqual(result.spectralEnergy,expected,'AbsTol',1e-12);
        end

        function laserDivergenceDeterminesWaist(testCase)
            [cfg,stack,folder]=singleLayerInput(testCase);
            stack.layers.singlelayer.method="reflectance";
            cfg.source.type="laser";cfg.source.laserSpecification="divergence";
            cfg.source.divergenceHalfAngleDeg=0.2;
            result=runCOEReflection(cfg,stack,folder);
            expected=cfg.wavelengths(2)/(pi*deg2rad(0.2));
            testCase.verifyEqual(result.illumination(2).beamRadius,expected,'RelTol',1e-12);
        end

        function outputRefinementRetainsAbsoluteEnergy(testCase)
            [cfg,stack,folder]=singleLayerInput(testCase);
            stack.layers.singlelayer.method="reflectance";
            coarse=runCOEReflection(cfg,stack,folder);
            cfg.outputSize=257;
            fine=runCOEReflection(cfg,stack,folder);
            testCase.verifyEqual(fine.metrics.Throughput,coarse.metrics.Throughput,'AbsTol',1e-12);
            testCase.verifyLessThan(max(abs(fine.metrics.CapturedThroughput-coarse.metrics.CapturedThroughput)),1e-3);
            testCase.verifySize(fine.kernel,[257,257,3]);
        end

        function gaussianBeamChangesShapeButNotMirrorReflectance(testCase)
            [cfg,stack,folder]=singleLayerInput(testCase);
            stack.layers.singlelayer.method="reflectance";
            stack.layers.singlelayer.reflectance=0.36;
            uniform=runCOEReflection(cfg,stack,folder);
            cfg.source.type="laser";
            cfg.source.waistRadius=40e-6;
            gaussian=runCOEReflection(cfg,stack,folder);
            testCase.verifyEqual(gaussian.metrics.Throughput,0.36*ones(3,1),'AbsTol',1e-12);
            testCase.verifyGreaterThan(max(abs(gaussian.kernel(:)-uniform.kernel(:))),1e-4);
            testCase.verifyEqual(gaussian.incidentPower,cfg.incidentIrradiance*gaussian.incidentArea,'AbsTol',1e-15);
            testCase.verifyLessThan(max(gaussian.metrics.ParsevalAbsoluteError),1e-12);
        end

        function finitePointSourceAveragesIntensitiesAndConservesPower(testCase)
            [cfg,stack,folder]=singleLayerInput(testCase);
            stack.layers.singlelayer.method="reflectance";
            stack.layers.singlelayer.reflectance=0.36;
            cfg.source.type="point";
            cfg.source.distance=0.1;
            cfg.source.diameter=0.2e-3;
            result=runCOEReflection(cfg,stack,folder);
            testCase.verifyEqual(result.metrics.Throughput,0.36*ones(3,1),'AbsTol',1e-12);
            testCase.verifyGreaterThan(result.illumination(1).modeCount,1);
            testCase.verifyLessThan(max(result.metrics.ParsevalAbsoluteError),1e-12);
            testCase.verifyEqual(mean(result.incidentIrradianceMap,'all'),cfg.incidentIrradiance,'AbsTol',1e-12);
        end

        function uniformIjpImageMatchesUniformFilm(testCase)
            [cfg,stack,folder]=singleLayerInput(testCase);
            stack.layerNames="IJP";
            stack.layers=struct('IJP',stack.layers.singlelayer);
            stack.layers.IJP.file="intentionally_missing.bmp";
            stack.layers.IJP.coverage="full";
            stack.layers.IJP.thickness=10e-6;
            result=runCOEReflection(cfg,stack,folder);
            expected=independentAirFilmReflectance(cfg.wavelengths,1.5,10e-6);
            testCase.verifyEqual(result.metrics.Throughput,expected(:),'AbsTol',1e-12);
            testCase.verifyTrue(all(result.layerGeometry.layerMasks.IJP,'all'));
            imwrite(uint8(255*ones(64)),fullfile(folder,'IJP.bmp'));
            stack.layers.IJP.coverage="image";stack.layers.IJP.file="IJP.bmp";
            fromImage=runCOEReflection(cfg,stack,folder);
            testCase.verifyEqual(fromImage.kernel,result.kernel,'AbsTol',1e-12);
            imwrite(false(64),fullfile(folder,'IJP.bmp'));
            stack.layers.IJP.invert=true;
            fromBlack=runCOEReflection(cfg,stack,folder);
            testCase.verifyTrue(all(fromBlack.layerGeometry.layerMasks.IJP,'all'));
            testCase.verifyEqual(fromBlack.config.coeStack.coefficientLookup, ...
                fromImage.config.coeStack.coefficientLookup,'AbsTol',1e-12);
            testCase.verifyEqual(fromBlack.spectralEnergy,fromImage.spectralEnergy,'AbsTol',1e-12);
        end

        function overlappingCfKeepsThicknessAndMatchesCharacteristicMatrix(testCase)
            [cfg,stack,folder]=overlappingCfInput(testCase);
            result=runCOEReflection(cfg,stack,folder);
            table=result.config.coeStack.regionTable;
            expected=overlappingCfReflectance(cfg.wavelengths,table,stack);
            testCase.verifyEqual(abs(result.config.coeStack.coefficientLookup).^2,expected,'AbsTol',1e-12);
            testCase.verifyEqual(result.metrics.Throughput, ...
                (table.AreaFraction.'*expected).','AbsTol',1e-12);
            testCase.verifyEqual(table.LocalHeightM, ...
                double(table.CF_R)*100e-9+double(table.CF_B)*130e-9+double(table.CF_G)*170e-9, ...
                'AbsTol',1e-18);
            testCase.verifyEqual(result.config.coeStack.referenceHeight,400e-9,'AbsTol',1e-18);
            testCase.verifyEqual(table.AreaFraction,0.25*ones(4,1),'AbsTol',1e-12);
        end

        function overlappingCfOrderIsCustomizable(testCase)
            [cfg,stack,folder]=overlappingCfInput(testCase);
            rbg=runCOEReflection(cfg,stack,folder);
            stack.layerNames=["CF_G","CF_B","CF_R"];
            gbr=runCOEReflection(cfg,stack,folder);
            table=gbr.config.coeStack.regionTable;
            expected=overlappingCfReflectance(cfg.wavelengths,table,stack);
            testCase.verifyEqual(abs(gbr.config.coeStack.coefficientLookup).^2,expected,'AbsTol',1e-12);
            testCase.verifyGreaterThan(max(abs(rbg.metrics.Throughput-gbr.metrics.Throughput)),1e-3);
        end

        function unifiedZeroHeightArrayUsesOnePropagation(testCase)
            [cfg,stack,folder]=singleLayerInput(testCase);
            cfg.heightMap=zeros(64);
            result=runCOEReflection(cfg,stack,folder);
            testCase.verifyFalse(result.heightPhaseEnabled);
            testCase.verifyEqual(result.propagationRuns,1);
            testCase.verifyEqual(result.kernel,result.baselineKernel);
        end

        function unifiedNonzeroHeightKeepsComparison(testCase)
            [cfg,stack,folder]=singleLayerInput(testCase);
            stack.layers.singlelayer.method="reflectance";
            cfg.heightMap=[zeros(64,32),ones(64,32)]*500e-9;
            result=runCOEReflection(cfg,stack,folder);
            testCase.verifyTrue(result.heightPhaseEnabled);
            testCase.verifyEqual(result.propagationRuns,2);
            testCase.verifyGreaterThan(max(abs(result.kernel(:)-result.baselineKernel(:))),1e-5);
        end

        function unifiedFreeStandingQuarterWaveFilm(testCase)
            [cfg,stack,folder]=singleLayerInput(testCase);
            stack.layers.singlelayer.thickness=550e-9/(4*1.5);
            result=runCOEReflection(cfg,stack,folder);
            expected=((1.5^2-1)/(1.5^2+1))^2;
            testCase.verifyEqual(result.metrics.Throughput(2),expected,'AbsTol',1e-12);
            testCase.verifyEqual(height(result.config.coeStack.materialTable),3);
            testCase.verifyEqual(result.config.coeStack.config.layerNames,"singlelayer");
        end

        function unifiedDirectReflectanceUsesAmplitudeAndPhase(testCase)
            [cfg,stack,folder]=singleLayerInput(testCase);
            stack.layers.singlelayer.method="reflectance";
            stack.layers.singlelayer.reflectance=[0.25,0.36,0.49];
            stack.layers.singlelayer.phaseDeg=[0,90,180];
            result=runCOEReflection(cfg,stack,folder);
            testCase.verifyEqual(result.metrics.Throughput,[0.25;0.36;0.49],'AbsTol',1e-12);
            testCase.verifyEqual(result.config.coeStack.coefficientLookup,[0.5,0.6i,-0.7],'AbsTol',1e-12);
        end

        function unifiedSameLevelOverlapIsRejected(testCase)
            [cfg,stack,folder]=singleLayerInput(testCase);
            stack.layers.other=stack.layers.singlelayer;
            stack.layerNames=["singlelayer","other"];
            stack.layerLevels=[1,1];
            testCase.verifyError(@() runCOEReflection(cfg,stack,folder),'COEReflection:SameLevelOverlap');
        end

        function unifiedFiniteFilmsRespectOrderAndAirExit(testCase)
            [cfg,stack,folder]=singleLayerInput(testCase);
            stack.layers.singlelayer.thickness=100e-9;
            stack.layers.other=stack.layers.singlelayer;
            stack.layers.other.n=2.1;
            stack.layers.other.thickness=130e-9;
            stack.layerNames=["singlelayer","other"];
            stack.layerLevels=[1,2];
            result=runCOEReflection(cfg,stack,folder);
            expected=independentAirFilmReflectance(cfg.wavelengths,[2.1,1.5],[130e-9,100e-9]);
            testCase.verifyEqual(result.metrics.Throughput,expected(:),'AbsTol',1e-12);
            stack.layerNames=["other","singlelayer"];
            reversed=runCOEReflection(cfg,stack,folder);
            testCase.verifyGreaterThan(max(abs(result.config.coeStack.coefficientLookup-reversed.config.coeStack.coefficientLookup)),0.01);
        end

        function bareMirrorAbsorberAndCommonReferencePlane(testCase)
            wavelengths = [450, 550, 650]*1e-9;
            mirror = coreReflectFilmStack(wavelengths, 1, zeros(0, 3), zeros(1, 0), -1, 0);
            absorber = coreReflectFilmStack(wavelengths, 1, zeros(0, 3), zeros(1, 0), 0, 0);
            testCase.verifyEqual(mirror, -ones(1, 3));
            testCase.verifyEqual(absorber, zeros(1, 3));
            airLayer = coreReflectFilmStack(wavelengths, 1, ones(1, 3), 100e-9, -1, 500e-9);
            bareAtReference = coreReflectFilmStack(wavelengths, 1, zeros(0, 3), zeros(1, 0), -1, 500e-9);
            testCase.verifyEqual(airLayer, bareAtReference, 'AbsTol', 2e-13);
        end

        function quarterWaveFilmCancelsKnownInterface(testCase)
            lambda = 550e-9;
            filmN = sqrt(1.5);
            depth = lambda/(4*filmN);
            coefficient = coreReflectFilmStack(lambda, 1, filmN, depth, -0.2, depth);
            testCase.verifyEqual(coefficient, 0, 'AbsTol', 2e-13);
            zeroDepth = coreReflectFilmStack(lambda, 1, filmN, 0, -0.2, 0);
            testCase.verifyEqual(zeroDepth, -0.2, 'AbsTol', 2e-13);
        end

        function losslessStackOnMirrorConservesUnitReflectance(testCase)
            wavelengths = (380:5:780)*1e-9;
            films = [1.7; 2.1].*ones(2, numel(wavelengths));
            depths = [120, 250]*1e-9;
            coefficient = coreReflectFilmStack(wavelengths, 1, films, depths, -1, sum(depths));
            testCase.verifyEqual(abs(coefficient).^2, ones(size(wavelengths)), 'AbsTol', 2e-13);
            lossy = coreReflectFilmStack(wavelengths, 1, 1.7 + 0.2i, 600e-9, -1, 600e-9);
            testCase.verifyLessThan(abs(lossy).^2, 1);
        end

        function opaqueFilmHidesUnderlyingTermination(testCase)
            wavelengths = [450, 550, 650]*1e-9;
            film = 1.8 + 5i;
            mirror = coreReflectFilmStack(wavelengths, 1, film, 1e-6, -1, 1e-6);
            absorber = coreReflectFilmStack(wavelengths, 1, film, 1e-6, 0, 1e-6);
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
            actual = coreReflectFilmStack(wavelengths, 1, indices, depths, 0, sum(depths));
            testCase.verifyEqual(actual, expected, 'AbsTol', 2e-13);
        end

        function bmOpeningDoesNotExpandAnodeReflectance(testCase)
            [cfg,stack,folder,data]=layerInput(testCase);
            stack.layers.IJP.thickness=0;
            stack.layers.BM.thickness=0;
            stack.layers.CF_R.thickness=0;
            stack.layers.CF_G.thickness=0;
            stack.layers.CF_B.thickness=0;
            result=runCOEReflection(cfg,stack,folder);
            testCase.verifyGreaterThan(nnz(~data.BM & ~data.Anode),0);
            testCase.verifyEqual(result.metrics.Throughput,mean(data.Anode,'all')*ones(3,1),'AbsTol',2e-13);
            model=result.config.coeStack;
            coefficient=reshape(model.coefficientLookup(double(model.regionMap(:))+1,1),64,64);
            testCase.verifyEqual(abs(coefficient),double(data.Anode),'AbsTol',2e-13);
            testCase.verifyEqual(result.kernel,result.baselineKernel,'AbsTol',2e-13);
        end

        function layerMATAndImagesRoundTripAndCSVNK(testCase)
            [cfg,stack,folder,data]=layerInput(testCase);
            expected=runCOEReflection(cfg,stack,folder);
            save(fullfile(folder,'layers.mat'),'-struct','data');
            stack.geometry.source="mat";
            fromMAT=runCOEReflection(cfg,stack,folder);
            testCase.verifyEqual(fromMAT.layerGeometry.layerMasks,expected.layerGeometry.layerMasks);
            testCase.verifyEqual(fromMAT.kernel,expected.kernel,'AbsTol',2e-13);
            stack.geometry.source="images";
            writematrix([450,1.9,0.35;550,1.8,0.35;650,1.7,0.35],fullfile(folder,'bm_nk.csv'));
            stack.layers.BM.nkFile="bm_nk.csv";
            fromImages=runCOEReflection(cfg,stack,folder);
            selected=fromImages.config.coeStack.materialTable.Material=="BM";
            testCase.verifyEqual(fromImages.config.coeStack.materialTable.n(selected),[1.9;1.8;1.7]);
            testCase.verifyEqual(fromImages.layerGeometry.layerMasks,data);
        end

        function rgbOverlapAndInvalidLayerOrderAreRejected(testCase)
            [cfg,stack,folder]=layerInput(testCase);
            stack.layers.CF_G.file="CF_R";
            testCase.verifyError(@() runCOEReflection(cfg,stack,folder),'COEReflection:SameLevelOverlap');
            stack.layers.CF_G.file="CF_G";
            stack.layerLevels=[1,2,3,4,3,4];
            testCase.verifyError(@() runCOEReflection(cfg,stack,folder),'COEReflection:LayerOrder');
        end

    end
end

function [repeated,complete]=repeatInputPair(testCase,inputPeriods)
    [cfg,stack,folder]=singleLayerInput(testCase);
    cfg.periods=[8,8]; cfg.outputSize=65;
    cfg.wavelengths=550e-9; cfg.spectralWeights=1; cfg.channelResponse=1;
    stack.layers.singlelayer.method="reflectance";
    stack.layers.singlelayer.reflectance=0.36;
    shape=inputPeriods([2,1])*cfg.samplesPerPeriod;
    [xx,yy]=meshgrid(0:shape(2)-1,0:shape(1)-1);
    unit=mod(floor(xx/3)+floor(yy/5),2)==0;
    imwrite(unit,fullfile(folder,'singlelayer.bmp'));
    stack.geometry.repeatToFit=true;
    repeated=runCOEReflection(cfg,stack,folder);
    imwrite(repmat(unit,cfg.periods([2,1])./inputPeriods([2,1])),fullfile(folder,'singlelayer.bmp'));
    stack.geometry.repeatToFit=false;
    complete=runCOEReflection(cfg,stack,folder);
end

function cfg = smallConfig()
    cfg = coeReflectionDefaults();
    cfg.periods = [4, 4]; cfg.samplesPerPeriod = 16; cfg.paddingFactor = 2; cfg.outputSize = 129;
    cfg.angleHalfRangeDeg = 3;
    cfg.wavelengths = [450, 550, 650]*1e-9; cfg.spectralWeights = [1, 1, 1]; cfg.channelResponse = eye(3);
end

function [cfg,stack,folder]=sourceConvolutionInput(testCase)
    [cfg,stack,folder]=singleLayerInput(testCase);
    cfg.wavelengths=550e-9; cfg.spectralWeights=1; cfg.channelResponse=1;
    cfg.outputSize=65; cfg.angleHalfRangeDeg=1.5; cfg.metricAngleDeg=1;
    stack.layers.singlelayer.method="reflectance";
    stack.layers.singlelayer.reflectance=0.36;
    cfg.source.type="point"; cfg.source.wavefront="paraxial";
    cfg.source.distance=0.05; cfg.source.samplesPerAxis=9;
    cfg.source.diameter=cfg.source.distance*deg2rad(3/64)*9;
end

function [cfg,stack,folder]=singleLayerInput(testCase)
    cfg=smallConfig();
    temporary=testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture());
    folder=string(temporary.Folder);
    imwrite(true(64),fullfile(folder,'singlelayer.bmp'));
    layer=struct('file',"singlelayer.bmp",'field',"singlelayer",'invert',false, ...
        'method',"nk",'n',1.5,'k',0,'thickness',0,'nkFile',"",'reflectance',1,'phaseDeg',180);
    stack=struct('incidentN',1,'exitN',1,'backgroundMode',"medium", ...
        'layerNames',"singlelayer",'layerLevels',1, ...
        'geometry',struct('source',"images",'matFile',"",'threshold',0.5), ...
        'layers',struct('singlelayer',layer));
end

function [cfg,stack,folder]=overlappingCfInput(testCase)
    [cfg,stack,folder]=singleLayerInput(testCase);
    stack.layerNames=["CF_R","CF_B","CF_G"];
    stack.layerLevels=[1,2,3];
    indices=[1.5,2.1,1.8];
    depths=[100,130,170]*1e-9;
    masks=false(64,64,3);
    masks(:,1:48,1)=true;
    masks(:,17:48,2)=true;
    masks(:,33:64,3)=true;
    for j=1:3
        name=stack.layerNames(j);
        stack.layers.(name)=stack.layers.singlelayer;
        stack.layers.(name).file=name+".png";
        stack.layers.(name).n=indices(j);
        stack.layers.(name).thickness=depths(j);
        imwrite(masks(:,:,j),fullfile(folder,name+".png"));
    end
    stack.layers.CF_B.k=0.05;
end

function expected=overlappingCfReflectance(wavelengths,regions,stack)
    expected=zeros(height(regions),numel(wavelengths));
    for row=1:height(regions)
        names=fliplr(stack.layerNames); % Independent characteristic matrices use top to bottom.
        present=logical(regions{row,cellstr(names)});
        layers=stack.layers;
        indices=arrayfun(@(name) layers.(name).n+1i*layers.(name).k,names(present));
        depths=arrayfun(@(name) layers.(name).thickness,names(present));
        expected(row,:)=independentAirFilmReflectance(wavelengths,indices,depths);
    end
end

function reflectance=independentAirFilmReflectance(wavelengths,indices,thickness)
    reflectance=zeros(size(wavelengths));
    for k=1:numel(wavelengths)
        matrix=eye(2);
        for j=1:numel(indices)
            phase=2*pi*indices(j)*thickness(j)/wavelengths(k);
            matrix=matrix*[cos(phase),-1i*sin(phase)/indices(j); ...
                -1i*indices(j)*sin(phase),cos(phase)];
        end
        admittance=(matrix(2,1)+matrix(2,2))/(matrix(1,1)+matrix(1,2));
        reflectance(k)=abs((1-admittance)/(1+admittance))^2;
    end
end

function expected=independentSourcePointSum(cfg,stack,folder)
    diameter=cfg.source.diameter;
    count=cfg.source.samplesPerAxis;
    axis=((0:count-1)+0.5)/count-0.5;
    [xx,yy]=meshgrid(axis*diameter);
    inside=hypot(xx,yy)<=diameter/2;
    xx=xx(inside); yy=yy(inside);
    expected=zeros(cfg.outputSize);total=0;
    cfg.source.diameter=0;
    for j=1:numel(xx)
        cfg.source.centerXY=[xx(j),yy(j)];
        point=runCOEReflection(cfg,stack,folder);
        power=1/point.illumination.normalizationScale^2;
        expected=expected+power*point.spectralEnergy;
        total=total+power;
    end
    expected=expected/total;
end

function [cfg,stack,folder,data]=layerInput(testCase)
    cfg=smallConfig();
    temporary=testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture());
    folder=string(temporary.Folder);
    names=["Anode","IJP","BM","CF_R","CF_G","CF_B"];
    base=struct('file',"",'field',"",'invert',false,'coverage',"image", ...
        'method',"nk",'n',1.5,'k',0,'thickness',100e-9,'nkFile',"",'reflectance',1,'phaseDeg',180);
    stack=struct('incidentN',1,'exitN',1,'backgroundMode',"absorber", ...
        'layerNames',names,'layerLevels',[1,2,3,4,4,4], ...
        'geometry',struct('source',"images",'matFile',"layers.mat",'threshold',0.5),'layers',struct());
    data.Anode=false(64);data.Anode(25:40,25:40)=true;
    data.IJP=true(64);
    data.BM=true(64);data.BM(17:48,17:48)=false;
    data.CF_R=false(64);data.CF_R(:,1:21)=true;
    data.CF_G=false(64);data.CF_G(:,22:43)=true;
    data.CF_B=false(64);data.CF_B(:,44:64)=true;
    for name=names
        stack.layers.(name)=base;
        stack.layers.(name).file=name;
        stack.layers.(name).field=name;
        imwrite(data.(name),fullfile(folder,name+".png"));
    end
    stack.layers.Anode.method="reflectance";stack.layers.Anode.thickness=0;
    stack.layers.IJP.thickness=10e-6;
    stack.layers.BM.k=0.35;
end
