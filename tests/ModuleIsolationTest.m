classdef ModuleIsolationTest < matlab.unittest.TestCase
    %ModuleIsolationTest Run copied feature folders without sibling dependencies.
    methods (Test)
        function coeCurrentD65PresetRunsFromCopiedModule(testCase)
            folder=portableModule(testCase,'topic_coe_reflection',false);
            script=fullfile(folder,'src','COEReflection_main.m');
            text=regexprep(fileread(script),'projectName = "[^"]*";', ...
                'projectName = "case03_PlaneWave_Pitch4";');
            file=fopen(script,'w','n','UTF-8'); fprintf(file,'%s',text); fclose(file);
            config=fullfile(folder,'Input','case03_PlaneWave_Pitch4','config.txt');
            text=strrep(fileread(config),'cfg.outputSize = 1801','cfg.outputSize = 257');
            file=fopen(config,'w','n','UTF-8'); fprintf(file,'%s',text); fclose(file);
            visibility=get(groot,'DefaultFigureVisible'); initial=findall(groot,'Type','figure');
            testCase.addTeardown(@()set(groot,'DefaultFigureVisible',visibility));
            testCase.addTeardown(@()closeNewFigures(initial)); set(groot,'DefaultFigureVisible','off');
            runScriptIsolated(script);
            output=fullfile(folder,'Output','case03_PlaneWave_Pitch4');
            saved=load(fullfile(output,'coe_reflection_results.mat'));
            testCase.verifyEqual(saved.cfg.wavelengths,(380:10:780)*1e-9,'AbsTol',1e-15);
            testCase.verifyEqual(saved.cfg.periods,[4,4]);
            testCase.verifyEqual(saved.cfg.propagationDistance,0.25);
            testCase.verifyFalse(saved.cfg.observer.enabled);
            testCase.verifySize(saved.colorData.withHeight.displaySRGB,[257,257,3]);
            testCase.verifyTrue(all(isfinite(saved.colorData.withHeight.XYZ),'all'));
            testCase.verifyEqual(saved.displayData.receiverDistanceView.distances,0.25);
            testCase.verifyTrue(isfile(fullfile(output,'coe_reflection_receiver_rgb_mm.png')));
        end

        function coeCoherentEyeRunsFromCopiedModule(testCase)
            folder=portableModule(testCase,'topic_coe_reflection');
            script=fullfile(folder,'src','COEReflection_main.m');
            text=regexprep(fileread(script),'projectName = "[^"]*";', ...
                'projectName = "case03_CoherentEye_Pitch4";');
            file=fopen(script,'w','n','UTF-8'); assert(file>=0); fprintf(file,'%s',text); fclose(file);
            config=fullfile(folder,'Input','case03_CoherentEye_Pitch4','config.txt');
            text=fileread(config);
            text=strrep(text,'cfg.outputSize = 2401','cfg.outputSize = 201');
            text=strrep(text,'cfg.observer.angleHalfRangeDeg = 4.5','cfg.observer.angleHalfRangeDeg = 0.15');
            text=strrep(text,'cfg.observer.pupilSamples = 1537','cfg.observer.pupilSamples = 129');
            text=strrep(text,'cfg.source.radialNodes = 7','cfg.source.radialNodes = 2');
            text=strrep(text,'cfg.source.azimuthNodes = 32','cfg.source.azimuthNodes = 8');
            file=fopen(config,'w','n','UTF-8'); assert(file>=0); fprintf(file,'%s',text); fclose(file);
            visibility=get(groot,'DefaultFigureVisible'); initial=findall(groot,'Type','figure');
            testCase.addTeardown(@() set(groot,'DefaultFigureVisible',visibility));
            testCase.addTeardown(@() closeNewFigures(initial)); set(groot,'DefaultFigureVisible','off');
            runScriptIsolated(script);
            output=fullfile(folder,'Output','case03_CoherentEye_Pitch4','CoherentEye');
            saved=load(fullfile(output,'coherent_eye_results.mat'));
            testCase.verifyEqual(saved.eyeImage.config.periods,[4,4]);
            testCase.verifyEqual(saved.eyeImage.focusDistance,0.5,'AbsTol',1e-14);
            testCase.verifyEqual(numel(saved.eyeImage.sourceWeights),16);
            testCase.verifyEqual(saved.eyeDisplay.mapping.mode,"legacy");
            testCase.verifyTrue(isfile(fullfile(output,'coherent_eye_comparison.png')));
            testCase.verifyGreaterThan(sum(saved.eyeImage.spectralEnergy,'all'),0);
        end

        function coeSmallSourceUsesBoundedHuePreservingDisplay(testCase)
            folder=portableModule(testCase,'topic_coe_reflection');
            script=fullfile(folder,'src','COEReflection_main.m');
            text=regexprep(fileread(script),'projectName = "[^"]*";', ...
                'projectName = "case03_FiniteEye_SmallSource";');
            file=fopen(script,'w','n','UTF-8'); assert(file>=0); fprintf(file,'%s',text); fclose(file);
            config=fullfile(folder,'Input','case03_FiniteEye_SmallSource','config.txt');
            text=strrep(fileread(config),'cfg.periods = [16, 16]','cfg.periods = [4, 4]');
            text=strrep(text,'cfg.outputSize = 1801','cfg.outputSize = 257');
            file=fopen(config,'w','n','UTF-8'); assert(file>=0); fprintf(file,'%s',text); fclose(file);
            visibility=get(groot,'DefaultFigureVisible'); initial=findall(groot,'Type','figure');
            testCase.addTeardown(@() set(groot,'DefaultFigureVisible',visibility));
            testCase.addTeardown(@() closeNewFigures(initial)); set(groot,'DefaultFigureVisible','off');
            runScriptIsolated(script);
            output=fullfile(folder,'Output','case03_FiniteEye_SmallSource','DistanceScaling');
            saved=load(fullfile(output,'distance_scaling_results.mat')); report=saved.report;
            view=report.views{report.baselineDisplayIndex};
            decoded=double(imread(fullfile(output,'observer_baseline_rgb.png')))/255;
            signal=view.channelSignal/view.binSolidAngle;
            unit=signal./max(max(signal,[],3),realmin);
            ratioError=cross(decoded,flipud(unit),3);

            testCase.verifyEqual(report.mapping.mode,"luminance");
            testCase.verifyEqual(view.sourceDiameter,0.0005,'AbsTol',1e-14);
            testCase.verifyEqual(view.config.periods,[4,4]);
            testCase.verifyLessThanOrEqual(max(decoded,[],'all'),0.7+0.5/255);
            testCase.verifyLessThanOrEqual(max(abs(ratioError),[],'all'),1/255+1e-10);
            testCase.verifyTrue(all(isfinite(view.spectralEnergy),'all'));
            testCase.verifyGreaterThan(sum(view.spectralEnergy,'all'),0);
        end

        function coeFiniteEyeRunsFromCopiedModule(testCase)
            folder=portableModule(testCase,'topic_coe_reflection');
            script=fullfile(folder,'src','COEReflection_main.m');
            text=regexprep(fileread(script),'projectName = "[^"]*";', ...
                'projectName = "case03_FiniteEye";');
            file=fopen(script,'w','n','UTF-8'); assert(file>=0); fprintf(file,'%s',text); fclose(file);
            config=fullfile(folder,'Input','case03_FiniteEye','config.txt');
            text=fileread(config);
            text=strrep(text,'cfg.periods = [80, 80]','cfg.periods = [4, 4]');
            text=strrep(text,'cfg.outputSize = 601','cfg.outputSize = 129');
            text=strrep(text,'[0.15, 0.20, 0.25, 0.30, 0.35]','[0.20, 0.25, 0.30]');
            file=fopen(config,'w','n','UTF-8'); assert(file>=0); fprintf(file,'%s',text); fclose(file);
            visibility=get(groot,'DefaultFigureVisible'); initial=findall(groot,'Type','figure');
            testCase.addTeardown(@() set(groot,'DefaultFigureVisible',visibility));
            testCase.addTeardown(@() closeNewFigures(initial)); set(groot,'DefaultFigureVisible','off');
            runScriptIsolated(script);
            output=fullfile(folder,'Output','case03_FiniteEye','DistanceScaling');
            saved=load(fullfile(output,'distance_scaling_results.mat'));

            testCase.verifyEqual(saved.report.baseConfig.source.diameter,5e-3,'AbsTol',1e-14);
            testCase.verifyEqual(saved.report.baseConfig.observer.pupilDiameter,3e-3,'AbsTol',1e-14);
            testCase.verifyEqual(saved.report.summary.GeometricDiameterMm([2,5]),[4;4],'AbsTol',1e-12);
            testCase.verifyTrue(all(isnan(saved.report.summary.MeasuredPeakDeg)));
            testCase.verifyLessThan(max(saved.report.summary.RayWeightInModel),1);
            testCase.verifyTrue(isfile(fullfile(output,'observer_participation_regions.png')));
            testCase.verifyTrue(isfile(fullfile(output,'observer_source_pupil_stages.png')));
        end

        function coeDistanceBenchmarkRunsFromCopiedModule(testCase)
            folder=prepareDistanceModule(testCase,true,false);
            runScriptIsolated(fullfile(folder,'src','COEReflection_main.m'));
            output=fullfile(folder,'Output','case03_DistanceScaling','DistanceScaling');
            saved=load(fullfile(output,'distance_scaling_results.mat'));

            testCase.verifyEqual(height(saved.report.summary),6);
            testCase.verifyEqual(saved.report.summary.AngularScale([2,5]),[0.5;0.5],'AbsTol',1e-14);
            testCase.verifyTrue(isfile(fullfile(output,'distance_scaling_comparison.png')));
            testCase.verifyTrue(isfile(fullfile(output,'distance_scaling_peaks.csv')));
        end

        function coeDistanceRgbViewRunsFromCopiedModule(testCase)
            folder=prepareDistanceModule(testCase,false,false);
            runScriptIsolated(fullfile(folder,'src','COEReflection_main.m'));
            output=fullfile(folder,'Output','case03_DistanceScaling');
            saved=load(fullfile(output,'coe_reflection_results.mat'),'result','displayData');

            testCase.verifyEqual(saved.result.distanceView.angularScale,0.5,'AbsTol',1e-14);
            testCase.verifySize(saved.displayData.distanceView.rgb,[129,129,3]);
            testCase.verifyTrue(isfile(fullfile(output,'coe_reflection_distance_view.png')));
            testCase.verifyTrue(isfile(fullfile(output,'coe_reflection_receiver_irradiance.png')));
        end

        function coeDistanceSpectralViewRunsFromCopiedModule(testCase)
            folder=prepareDistanceModule(testCase,false,true);
            runScriptIsolated(fullfile(folder,'src','COEReflection_main.m'));
            output=fullfile(folder,'Output','case03_DistanceScaling');
            saved=load(fullfile(output,'coe_reflection_results.mat'),'result','displayData');

            testCase.verifySize(saved.result.distanceView.spectralEnergy,[129,129,2]);
            testCase.verifySize(saved.displayData.distanceView.rgb,[129,129,3]);
            testCase.verifyTrue(isfile(fullfile(output,'coe_reflection_distance_view.png')));
        end

        function coeFourPatternsRunFromCopiedModule(testCase)
            folder=portableModule(testCase,'topic_coe_reflection');
            scriptPath=fullfile(folder,'src','COEReflection_main.m');
            scriptText=strrep(fileread(scriptPath),'runFourPatternComparison = false;', ...
                'runFourPatternComparison = true;');
            scriptText=regexprep(scriptText,'projectName = "[^"]*";', ...
                'projectName = "case01";');
            scriptFile=fopen(scriptPath,'w','n','UTF-8'); assert(scriptFile>=0);
            fprintf(scriptFile,'%s',scriptText); fclose(scriptFile);
            index=floor(((0:31)+0.5)*500/32)+1;
            for name=["case01","case02","case03","case04"]
                projectDir=fullfile(folder,'Input',name);
                for stem=["anode","BM"]
                    imageData=imread(fullfile(projectDir,stem+".png"));
                    imwrite(imageData(index,index,:),fullfile(projectDir,stem+".png"));
                end
                file=fullfile(projectDir,'config.txt'); text=fileread(file);
                text=strrep(text,'cfg.periods = [8, 8]','cfg.periods = [4, 4]');
                text=strrep(text,'cfg.samplesPerPeriod = 250','cfg.samplesPerPeriod = 16');
                text=strrep(text,'stack.geometry.inputSamplesPerPeriod = 250', ...
                    'stack.geometry.inputSamplesPerPeriod = 16');
                text=strrep(text,'cfg.outputSize = 1201','cfg.outputSize = 129');
                text=strrep(text,'cfg.angleHalfRangeDeg = 6','cfg.angleHalfRangeDeg = 3');
                handle=fopen(file,'w','n','UTF-8'); assert(handle>=0);
                fprintf(handle,'%s',text); fclose(handle);
            end
            visibility=get(groot,'DefaultFigureVisible');
            testCase.addTeardown(@() set(groot,'DefaultFigureVisible',visibility));
            initialFigures=findall(groot,'Type','figure');
            testCase.addTeardown(@() closeNewFigures(initialFigures));
            set(groot,'DefaultFigureVisible','off');
            runScriptIsolated(fullfile(folder,'src','COEReflection_main.m'));
            output=fullfile(folder,'Output','FourPatterns');
            saved=load(fullfile(output,'four_patterns_results.mat')); report=saved.report;
            testCase.verifyEqual(report.caseNames,["case01","case02","case03","case04"]);
            testCase.verifyEqual(height(report.summary),24);
            testCase.verifyTrue(isfile(fullfile(output,'four_patterns_comparison.png')));
            for j=1:4
                testCase.verifyEqual(report.local{j}.cfg.periods,[2,2]);
                testCase.verifyEqual(report.array{j}.repeatCounts,[2,2]);
                testCase.verifyEqual(mean(report.patterns{j},'all'), ...
                    report.local{j}.metrics.Throughput(2),'AbsTol',1e-12);
                weights=reshape(report.local{j}.cfg.spectralWeights,1,1,3);
                signal=report.local{j}.spectralEnergy.*weights;
                expected=(log1p(report.mapping.alpha*min(signal/report.mapping.referenceSignal,1)) ...
                    /log1p(report.mapping.alpha)).^report.mapping.gamma;
                file=fullfile(folder,'Output',report.caseNames(j),'local_rgb.png');
                decoded=double(imread(file))/255;
                testCase.verifyEqual(decoded,flipud(expected),'AbsTol',0.5/255+eps);
            end
        end

        function coePlaneWaveComparisonRunsFromCopiedModule(testCase)
            folder=portableModule(testCase,'topic_coe_reflection');
            scriptPath=fullfile(folder,'src','COEReflection_main.m');
            scriptText=regexprep(fileread(scriptPath),'projectName = "[^"]*";', ...
                'projectName = "F192_AnodeBM_PlaneRGB";');
            scriptFile=fopen(scriptPath,'w','n','UTF-8');
            assert(scriptFile>=0); fprintf(scriptFile,'%s',scriptText); fclose(scriptFile);
            projectDir=fullfile(folder,'Input','F192_AnodeBM_PlaneRGB');
            % 临时夹具使用64×64原生图片；生产真实图样不变。
            index=floor(((0:63)+0.5)*1000/64)+1;
            for stem=["anode","BM"]
                imageData=imread(fullfile(projectDir,stem+".png"));
                imwrite(imageData(index,index,:),fullfile(projectDir,stem+".png"));
            end
            file=fullfile(projectDir,'config.txt');
            text=fileread(file);
            text=strrep(text,'cfg.periods = [16, 16]','cfg.periods = [8, 8]');
            text=strrep(text,'cfg.samplesPerPeriod = 250','cfg.samplesPerPeriod = 16');
            text=strrep(text,'stack.geometry.inputSamplesPerPeriod = 250', ...
                'stack.geometry.inputSamplesPerPeriod = 16');
            text=strrep(text,'cfg.outputSize = 1201','cfg.outputSize = 129');
            text=strrep(text,'cfg.angleHalfRangeDeg = 6','cfg.angleHalfRangeDeg = 3');
            handle=fopen(file,'w','n','UTF-8');
            assert(handle>=0); fprintf(handle,'%s',text); fclose(handle);
            visibility=get(groot,'DefaultFigureVisible');
            testCase.addTeardown(@() set(groot,'DefaultFigureVisible',visibility));
            initialFigures=findall(groot,'Type','figure');
            testCase.addTeardown(@() closeNewFigures(initialFigures));
            set(groot,'DefaultFigureVisible','off');
            runScriptIsolated(fullfile(folder,'src','COEReflection_main.m'));
            output=fullfile(folder,'Output','F192_AnodeBM_PlaneRGB','PlaneWaveComparison');
            saved=load(fullfile(output,'plane_wave_comparison.mat'));
            report=saved.report;
            testCase.verifyEqual(report.local.cfg.periods,[4,4]);
            testCase.verifyEqual(report.array.cfg.periods,[8,8]);
            testCase.verifyEqual(report.array.repeatCounts,[2,2]);
            testCase.verifyEqual(report.local.coordinate,report.array.coordinate);
            testCase.verifyEqual(report.local.metrics.Throughput, ...
                report.array.metrics.Throughput,'AbsTol',1e-12);
            testCase.verifyEqual(report.array.incidentPower/report.local.incidentPower,4,'AbsTol',1e-12);
            testCase.verifyEqual(mean(abs(report.local.effectiveReflection).^2,'all'), ...
                report.local.metrics.Throughput(2),'AbsTol',1e-12);
            expected=(log1p(report.mapping.alpha*min( ...
                report.array.channelSignal/report.mapping.referenceSignal,1)) ...
                /log1p(report.mapping.alpha)).^report.mapping.gamma;
            testCase.verifyEqual(report.array.rgb,expected,'AbsTol',1e-14);
            testCase.verifySize(report.local.rgb,[129,129,3]);
            testCase.verifyTrue(isfile(fullfile(output,'plane_wave_comparison.png')));
            testCase.verifyTrue(isfile(fullfile(output,'geometry_and_profile.png')));
        end

        function coeSourceConvolutionRunsFromCopiedModule(testCase)
            folder=portableModule(testCase,'topic_coe_reflection');
            cfg=portableConfig(coeReflectionDefaults());
            cfg.complexCoefficient=0.5;
            cfg.source.type="point";
            cfg.source.wavefront="plane";
            cfg.source.integrationMethod="angularConvolution";
            cfg.source.distance=0.05; cfg.source.diameter=0.2e-3;
            value=runCOEReflection(cfg);
            testCase.verifyEqual(value.metrics.Throughput,0.25*ones(3,1),'AbsTol',1e-12);
            testCase.verifyEqual(value.sourceConvolution(1).kernelMass,1,'AbsTol',1e-14);
            testCase.verifyEqual(which('runCOEReflection'),fullfile(folder,'src','runCOEReflection.m'));
        end

        function coeMultilayerMainRunsWithOnlyCopiedModule(testCase)
            folder = portableModule(testCase, 'topic_coe_reflection');
            scriptPath = fullfile(folder, 'src', 'COEReflection_main.m');
            prepareSmallScript(testCase, scriptPath);

            runScriptIsolated(scriptPath);
            saved = load(fullfile(folder, 'Output', 'coe_reflection_results.mat'));
            testCase.verifyEqual(saved.stack.layerLevels, [1,2,3,4,5,6]);
            testCase.verifyEqual(saved.stack.layerNames, ["Anode","IJP","BM","CF_R","CF_B","CF_G"]);
            testCase.verifyEqual(saved.cfg.heightMap, 0);
            testCase.verifyEqual(saved.result.kernel, saved.baseline.kernel, 'AbsTol', 2e-13);
            testCase.verifyGreaterThan(nnz(saved.stackGeometry.layerMasks.CF_R & ~saved.stackGeometry.layerMasks.Anode & ~saved.stackGeometry.layerMasks.BM), 0);
            testCase.verifyTrue(isfile(fullfile(folder, 'Output', 'coe_reflection_layer_patterns.mat')));
            testCase.verifyTrue(isfile(fullfile(folder, 'Output', 'coe_reflection_stack_materials.csv')));
            expected = saved.cfg.coeStack.regionTable.AreaFraction.' ...
                *abs(saved.cfg.coeStack.coefficientLookup(saved.cfg.coeStack.regionTable.RegionCode + 1, :)).^2;
            testCase.verifyEqual(saved.result.metrics.Throughput, expected.', 'AbsTol', 2e-13);
        end

        function coeD65ScriptProducesColorAndAngularResults(testCase)
            folder = portableModule(testCase, 'topic_coe_reflection');
            scriptPath = fullfile(folder, 'src', 'COEReflection_main.m');
            prepareSmallScript(testCase, scriptPath, "D65");

            runScriptIsolated(scriptPath);

            outputDir = fullfile(folder, 'Output');
            saved = load(fullfile(outputDir, 'coe_reflection_results.mat'));
            testCase.verifyEqual(saved.cfg.wavelengths, (380:5:780)*1e-9, 'AbsTol', 1e-21);
            testCase.verifyEqual(saved.spectrum.whiteXY, [0.3127, 0.3290], 'AbsTol', 1e-4);
            testCase.verifyEqual(saved.displayData.withHeightRGB, saved.colorData.withHeight.displaySRGB);
            testCase.verifySize(saved.colorData.withHeight.XYZ, [129, 129, 3]);
            testCase.verifyTrue(isfile(fullfile(outputDir, 'coe_reflection_diffraction_spots.csv')));
            testCase.verifyNotEmpty(saved.spotMetrics);
            testCase.verifyEqual(unique(saved.spotMetrics.StateIndex),1);
            testCase.verifyFalse(saved.heightPhaseEnabled);
            testCase.verifyEqual(saved.result.propagationRuns,1);
            testCase.verifyTrue(isfile(fullfile(outputDir, 'coe_reflection_cumulative_energy.png')));
            testCase.verifyTrue(isfile(fullfile(outputDir, 'coe_reflection_spectrum.csv')));
            obs = readtable(fullfile(outputDir, 'coe_reflection_observation_metrics.csv'));
            testCase.verifyEqual(obs.ViewRadiantReflectancePercent(1), ...
                100*saved.withAngular.radiantViewReflectance, 'AbsTol', 1e-10);
            testCase.verifyEqual(height(obs),1);
            testCase.verifyEqual(obs.TotalRadiantReflectancePercent,100*saved.withAngular.radiantTotalReflectance,'AbsTol',1e-10);
            cumulative = readtable(fullfile(outputDir, 'coe_reflection_cumulative_energy.csv'));
            testCase.verifyGreaterThanOrEqual(diff(cumulative.RadiantPercent), -1e-12);
            testCase.verifyEqual(cumulative.RadiantPercent(cumulative.RadiusDeg == 2.5), ...
                obs.ViewRadiantReflectancePercent(1), 'AbsTol', 1e-10);
        end

        function coeMainScriptProducesResults(testCase)
            folder = portableModule(testCase, 'topic_coe_reflection');
            scriptPath = fullfile(folder, 'src', 'COEReflection_main.m');
            prepareSmallScript(testCase, scriptPath, "RGB", "height_test.mat");
            % 用户输入非零高度 MAT；不再依赖生产代码生成条纹。
            heightValue=687.5e-9;
            heightMap=repmat([1,1,0,0],64,16)*heightValue;
            save(fullfile(folder,'Input','height_test.mat'),'heightMap');

            runScriptIsolated(scriptPath);

            outputDir = fullfile(folder, 'Output');
            saved = load(fullfile(outputDir, 'coe_reflection_results.mat'));
            expected = runCOEReflection(saved.cfg);
            testCase.verifyEqual(saved.result.kernel, expected.kernel, 'AbsTol', 2e-13);
            testCase.verifyEqual(saved.result.metrics, expected.metrics);
            testCase.verifyEqual(unique(saved.cfg.heightMap), [0; 687.5e-9]);
            testCase.verifyEqual(mean(saved.cfg.heightMap(:) > 0), 0.5);
            testCase.verifyEqual(saved.cfg.heightMap(1, 1:8), ...
                [1, 1, 0, 0, 1, 1, 0, 0]*687.5e-9);
            testCase.verifyGreaterThan(max(abs(saved.result.kernel(:) - saved.baseline.kernel(:))), 1e-6);
            testCase.verifyEqual(saved.result.metrics.Throughput, saved.baseline.metrics.Throughput, ...
                'AbsTol', 2e-13);
            testCase.verifyLessThan(saved.result.metrics.ZeroOrderRetention(2), 0.1);
            metrics = readtable(fullfile(outputDir, 'coe_reflection_metrics.csv'));
            testCase.verifyEqual(metrics.ZeroOrderRetention, saved.result.metrics.ZeroOrderRetention, 'AbsTol', 1e-10);
            testCase.verifyTrue(isfile(fullfile(outputDir, 'coe_reflection_baseline_metrics.csv')));
            comparison = readtable(fullfile(outputDir, 'coe_reflection_comparison_metrics.csv'), ...
                'TextType', 'string', 'Encoding', 'UTF-8');
            testCase.verifyEqual(comparison.ConeHalfAngleDeg, 2.5*ones(6, 1));
            testCase.verifyEqual(comparison.ConeFullAngleDeg, 5*ones(6, 1));
            testCase.verifyEqual(comparison.ReflectionHazeProxyPercent, ...
                100*[saved.baseline.metrics.OutsideAngleFraction; saved.result.metrics.OutsideAngleFraction], ...
                'AbsTol', 1e-10);
            testCase.verifyEqual(comparison.HeightPhaseState, ...
                [repmat("不叠加高度相位", 3, 1); repmat("叠加高度相位", 3, 1)]);
            testCase.verifyGreaterThan(max(abs(saved.displayData.withoutHeightRGB(:) ...
                - saved.displayData.withHeightRGB(:))), 0.01);
            testCase.verifyGreaterThanOrEqual(min(saved.displayData.withHeightRGB(:)), 0);
            testCase.verifyLessThanOrEqual(max(saved.displayData.withHeightRGB(:)), 1);
            for state = ["without", "with"]
                info = imfinfo(fullfile(outputDir, "coe_reflection_rgb_" + state + "_height_phase.png"));
                testCase.verifyEqual(info.ColorType, 'truecolor');
            end
            for suffix = ["geometry", "diffraction", "metrics"]
                info = imfinfo(fullfile(outputDir, "coe_reflection_" + suffix + ".png"));
                testCase.verifyGreaterThan(info.Width, 0);
                testCase.verifyGreaterThan(info.Height, 0);
            end
        end

        function coeMainScriptReadsDesignInput(testCase)
            folder = portableModule(testCase, 'topic_coe_reflection');
            scriptPath = fullfile(folder, 'src', 'COEReflection_main.m');
            prepareSmallScript(testCase, scriptPath, "RGB", "design.mat");
            heightMap = [zeros(64, 32), ones(64, 32)]*550e-9/(4*0.20);
            complexCoefficient = 0.5;
            save(fullfile(folder, 'Input', 'design.mat'), 'heightMap', 'complexCoefficient');

            runScriptIsolated(scriptPath);

            saved = load(fullfile(folder, 'Output', 'coe_reflection_results.mat'));
            testCase.verifyEqual(saved.cfg.mask, double(saved.stackGeometry.unionMask));
            testCase.verifyEqual(saved.cfg.heightMap, heightMap);
            testCase.verifyEqual(saved.cfg.complexCoefficient, complexCoefficient);
            testCase.verifyEqual(saved.baseline.config.heightMap, 0);
            expected = runCOEReflection(saved.cfg);
            testCase.verifyEqual(saved.result.kernel, expected.kernel, 'AbsTol', 2e-13);
            testCase.verifyEqual(saved.result.metrics, expected.metrics);
            testCase.verifyLessThan(saved.result.metrics.ZeroOrderRetention(2), 1);
            % RGB 映射必须在两种状态间共享参考值，不能分别归一化掩盖差异。
            mapping = saved.displayData.rgbMapping;
            expectedRGB = (log1p(mapping.alpha*min(saved.result.channelSignal/mapping.referenceSignal, 1)) ...
                /log1p(mapping.alpha)).^mapping.gamma;
            expectedWithout = (log1p(mapping.alpha*min(saved.baseline.channelSignal/mapping.referenceSignal, 1)) ...
                /log1p(mapping.alpha)).^mapping.gamma;
            testCase.verifyEqual(saved.displayData.withHeightRGB, expectedRGB, 'AbsTol', 2e-13);
            testCase.verifyEqual(saved.displayData.withoutHeightRGB, expectedWithout, 'AbsTol', 2e-13);
        end

        function coeHeightDefaultsToZero(testCase)
            folder = portableModule(testCase, 'topic_coe_reflection');
            scriptPath = fullfile(folder, 'src', 'COEReflection_main.m');
            prepareSmallScript(testCase, scriptPath);

            runScriptIsolated(scriptPath);

            saved = load(fullfile(folder, 'Output', 'coe_reflection_results.mat'));
            testCase.verifyFalse(saved.heightPhaseEnabled);
            testCase.verifyEqual(saved.result.propagationRuns,1);
            testCase.verifyEqual(height(saved.observationMetrics),1);
            testCase.verifyEqual(height(saved.comparisonMetrics),3);
            testCase.verifyFalse(isfile(fullfile(folder,'Output','coe_reflection_baseline_metrics.csv')));
            testCase.verifyTrue(isfile(fullfile(folder,'Output','coe_reflection_rgb.png')));
            farField=findall(groot,'Type','figure','Name','COE 远场反射');
            testCase.verifyEqual(numel(findall(farField,'Type','axes')),4);
            testCase.verifyEqual(saved.cfg.heightMap, 0);
            testCase.verifyEqual(saved.result.kernel, saved.baseline.kernel, 'AbsTol', 2e-13);
            testCase.verifyEqual(saved.result.metrics.ZeroOrderRetention, ones(3, 1), 'AbsTol', 2e-13);
        end

        function coeFolderRunsFromTemporaryLocation(testCase)
            folder = portableModule(testCase, 'topic_coe_reflection');
            cfg = portableConfig(coeReflectionDefaults());
            cfg.complexCoefficient = 0.5;

            initialPath = path;
            result = runCOEReflection(cfg);
            testCase.verifyEqual(path, initialPath);

            testCase.verifyEqual(result.metrics.Throughput, 0.25*ones(3, 1), 'AbsTol', 2e-13);
            testCase.verifyEqual(which('runCOEReflection'), fullfile(folder, 'src', 'runCOEReflection.m'));
        end

        function udcFolderCarriesItsDefaultImage(testCase)
            folder = portableModule(testCase, 'topic_udc_imaging');
            cfg = portableConfig(udcImagingDefaults());
            cfg.complexCoefficient = 0.5;
            imageInfo = imfinfo(fullfile(folder, 'Input', 'USAF500.png'));

            initialPath = path;
            [result, images] = runUDCImaging(cfg);
            testCase.verifyEqual(path, initialPath);

            testCase.verifyEqual(result.channelThroughput, [0.25, 0.25, 0.25], 'AbsTol', 2e-13);
            testCase.verifySize(images.inputLinear, [imageInfo.Height, imageInfo.Width, 3]);
            testCase.verifyEqual(which('runUDCImaging'), fullfile(folder, 'src', 'runUDCImaging.m'));
        end

        function apertureScriptRunsWithOnlyItsFolder(testCase)
            folder = portableModule(testCase, 'topic_aperture_fraunhofer');
            scriptPath = fullfile(folder, 'src', 'CustomAperture_fraunhofer_fft.m');
            prepareSmallScript(testCase, scriptPath);

            runScriptIsolated(scriptPath);

            testCase.verifyTrue(isfile(fullfile(folder, 'Output', 'four_aperture_rgb_preview.png')));
            testCase.verifyTrue(isfile(fullfile(folder, 'Output', 'four_aperture_gray_preview.png')));
        end

        function aperturePreviewRunsWithOnlyItsSourceFolder(testCase)
            folder = portableModule(testCase, 'topic_aperture_fraunhofer');
            scriptPath = fullfile(folder, 'src', 'CustomAperture_four_aperture_preview.m');
            prepareSmallScript(testCase, scriptPath);

            runScriptIsolated(scriptPath);

            testCase.verifyTrue(isfile(fullfile(folder, 'Output', 'four_aperture_preview_data.mat')));
        end

        function phaseScriptRunsWithOnlyItsFolder(testCase)
            folder = portableModule(testCase, 'topic_phase_layer');
            scriptPath = fullfile(folder, 'src', 'CustomApertureWithPhase_fraunhofer_fft.m');
            prepareSmallScript(testCase, scriptPath);

            runScriptIsolated(scriptPath);

            testCase.verifyTrue(isfile(fullfile(folder, 'Output', 'four_aperture_with_phase_rgb_preview.png')));
            testCase.verifyTrue(isfile(fullfile(folder, 'Output', 'aperiodic_phase_map_2d_3d.png')));
        end

        function continuousApertureUsesBundledSampler(testCase)
            folder = portableModule(testCase, 'topic_aperture_fraunhofer');
            scriptPath = fullfile(folder, 'src', 'CustomAperture_fraunhofer_fft.m');
            prepareSmallScript(testCase, scriptPath);
            selectContinuousAperture(scriptPath);

            runScriptIsolated(scriptPath);

            testCase.verifyTrue(isfile(fullfile(folder, 'Output', 'four_aperture_rgb_preview.png')));
        end

        function pixelGuiStartsFromTemporaryLocation(testCase)
            folder = portableModule(testCase, 'topic_pixel_unfolding');
            visibility = get(groot, 'DefaultFigureVisible');
            testCase.addTeardown(@() set(groot, 'DefaultFigureVisible', visibility));
            initialFigures = findall(groot, 'Type', 'figure');
            testCase.addTeardown(@() closeNewFigures(initialFigures));
            set(groot, 'DefaultFigureVisible', 'off');

            OLED_PixelUnfolding_Tian_Batch_GUI();

            testCase.verifyEqual(which('OLED_PixelUnfolding_Tian_Batch_GUI'), ...
                fullfile(folder, 'src', 'OLED_PixelUnfolding_Tian_Batch_GUI.m'));
            testCase.verifyNotEmpty(setdiff(findall(groot, 'Type', 'figure'), initialFigures));
        end

        function phaseLoaderReadsLocalOutput(testCase)
            folder = portableModule(testCase, 'topic_phase_layer');
            resultPath = writeMinimalDoeResult(folder);

            design = loadBestMultiRadiusPhaseDesign();

            testCase.verifyEqual(design.sourcePath, resultPath);
            testCase.verifyEqual(design.phaseIslandRadiusList, [1, 1.5]*1e-6, 'AbsTol', 1e-18);
            testCase.verifyEqual(design.phaseHeight, 220e-9, 'AbsTol', 1e-18);
        end
    end
end

function folder=prepareDistanceModule(testCase,benchmark,spectral)
    folder=portableModule(testCase,'topic_coe_reflection');
    scriptPath=fullfile(folder,'src','COEReflection_main.m');
    text=regexprep(fileread(scriptPath),'projectName = "[^"]*";', ...
        'projectName = "case03_DistanceScaling";');
    file=fopen(scriptPath,'w','n','UTF-8'); assert(file>=0);
    fprintf(file,'%s',text); fclose(file);
    project=fullfile(folder,'Input','case03_DistanceScaling');
    index=floor(((0:31)+0.5)*500/32)+1;
    for name=["anode","BM"]
        image=imread(fullfile(project,name+".png"));
        imwrite(image(index,index,:),fullfile(project,name+".png"));
    end
    text=fileread(fullfile(project,'config.txt'));
    text=strrep(text,'cfg.samplesPerPeriod = 250','cfg.samplesPerPeriod = 16');
    text=strrep(text,'stack.geometry.inputSamplesPerPeriod = 250', ...
        'stack.geometry.inputSamplesPerPeriod = 16');
    text=strrep(text,'cfg.outputSize = 601','cfg.outputSize = 129');
    text=strrep(text,'cfg.angleHalfRangeDeg = 6','cfg.angleHalfRangeDeg = 3');
    text=strrep(text,'[0.15, 0.20, 0.25, 0.30, 0.35]','[0.20, 0.25, 0.30]');
    if ~benchmark, text=strrep(text,'cfg.observer.runBenchmark = true','cfg.observer.runBenchmark = false'); end
    if spectral
        text=strrep(text,'spectrumMode = "RGB"','spectrumMode = "D65"');
        text=[text,newline,'wavelengthRangeNm = [540, 560]',newline,'wavelengthStepNm = 20',newline];
    end
    file=fopen(fullfile(project,'config.txt'),'w','n','UTF-8'); assert(file>=0);
    fprintf(file,'%s',text); fclose(file);
    visibility=get(groot,'DefaultFigureVisible');
    testCase.addTeardown(@() set(groot,'DefaultFigureVisible',visibility));
    initialFigures=findall(groot,'Type','figure');
    testCase.addTeardown(@() closeNewFigures(initialFigures));
    set(groot,'DefaultFigureVisible','off');
end

function folder = portableModule(testCase, moduleName, legacyConfig)
    if nargin<3, legacyConfig=true; end
    root = fileparts(fileparts(mfilename('fullpath')));
    temporary = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture());
    folder = fullfile(temporary.Folder, moduleName);
    mkdir(folder);
    % Never copy historical outputs: output checks must observe newly generated files.
    for part = ["Input", "src", "docs", "dist"]
        copyfile(fullfile(root, moduleName, part), fullfile(folder, part));
    end
    % Legacy branch tests use frozen parameters, independent of current user presets.
    if moduleName=="topic_coe_reflection" && legacyConfig
        fixtures=dir(fullfile(root,'tests','fixtures','coe_legacy_configs','*.legacy.txt'));
        for j=1:numel(fixtures)
            project=erase(string(fixtures(j).name),'.legacy.txt');
            copyfile(fullfile(fixtures(j).folder,fixtures(j).name),fullfile(folder,'Input',project,'config.txt'));
        end
    end
    mkdir(fullfile(folder, 'Output'));
    testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(folder, 'src')));
end

function cfg = portableConfig(cfg)
    cfg.periods = [4, 4];
    cfg.samplesPerPeriod = 16;
    cfg.paddingFactor = 2;
    cfg.pupilRadius = Inf;
    cfg.outputSize = 129;
    cfg.sensorPitch = 10e-6;
    cfg.angleHalfRangeDeg = 3;
end

function prepareSmallScript(testCase, scriptPath, spectrumMode, inputFile)
    arguments
        testCase
        scriptPath
        spectrumMode (1,1) string = "RGB"
        inputFile (1,1) string = ""
    end
    % Only the temporary copy is changed; production simulation defaults stay intact.
    text = fileread(scriptPath);
    if endsWith(scriptPath, 'COEReflection_main.m')
        text = regexprep(text, '(?m)^projectName = "[^"]*";', 'projectName = "";');
        text = regexprep(text, '(?m)^useProjectConfig = [^;\r\n]+;', 'useProjectConfig = false;');
        text = regexprep(text, '(?m)^inputFile = "[^"]*";', sprintf('inputFile = "%s";', inputFile));
        text = regexprep(text, '(?m)^spectrumMode = "[^"]*";', sprintf('spectrumMode = "%s";', spectrumMode));
        text = regexprep(text, '(?m)^wavelengthStepNm = [^;\r\n]+;', 'wavelengthStepNm = 5;');
        text = regexprep(text, '(?m)^cfg\.periods = [^;\r\n]+;', 'cfg.periods = [4, 4];');
        text = regexprep(text, '(?m)^cfg\.samplesPerPeriod = [^;\r\n]+;', 'cfg.samplesPerPeriod = 16;');
        text = regexprep(text, '(?m)^cfg\.outputSize = [^;\r\n]+;', 'cfg.outputSize = 129;');
        % At 50 um pitch, 16 samples and 380 nm, Nyquist is about 3.48 deg.
        % The outer pixel edge of this 129-pixel grid is about 3.02 deg.
        text = regexprep(text, '(?m)^cfg\.angleHalfRangeDeg = [^;\r\n]+;', 'cfg.angleHalfRangeDeg = 3;');
        inputFolder = fullfile(fileparts(fileparts(scriptPath)), 'Input');
        for stem = ["anode","IJP","BM","CF_R","CF_G","CF_B"]
            source = imread(fullfile(inputFolder,'Case1',stem+".png"));
            imwrite(source(1:4:256,1:4:256),fullfile(inputFolder,stem+".png"));
        end
    else
        text = strrep(text, 'clear, clc, close all', 'clear');
        text = strrep(text, 'nn = 1000;', 'nn = 24;');
        text = strrep(text, 'nn = 700;', 'nn = 24;');
        text = strrep(text, 'useBestDoeDesign = true;', 'useBestDoeDesign = false;');
        text = strrep(text, 'phaseFillFactor = 0.40;', 'phaseFillFactor = 0.01;');
    end
    file = fopen(scriptPath, 'w', 'n', 'UTF-8');
    assert(file >= 0, 'Cannot write temporary script.');
    closeFile = onCleanup(@() fclose(file));
    fprintf(file, '%s', text);
    originalVisibility = get(groot, 'DefaultFigureVisible');
    testCase.addTeardown(@() set(groot, 'DefaultFigureVisible', originalVisibility));
    set(groot, 'DefaultFigureVisible', 'off');
    initialFigures = findall(groot, 'Type', 'figure');
    testCase.addTeardown(@() closeNewFigures(initialFigures));
    initialRng = rng;
    testCase.addTeardown(@() rng(initialRng));
end

function runScriptIsolated(scriptPath)
    % Legacy scripts clear their caller workspace; keep that away from TestCase.
    run(scriptPath);
end

function closeNewFigures(initialFigures)
    created = setdiff(findall(groot, 'Type', 'figure'), initialFigures);
    for figureHandle = created(:).'
        close(figureHandle);
    end
end

function resultPath = writeMinimalDoeResult(folder)
    outputDir = fullfile(folder, 'Output');
    if ~isfolder(outputDir)
        mkdir(outputDir);
    end
    results = table("[1 1.5]", 0.1, 1.0, 5, 0.2, 220, ...
        'VariableNames', {'RadiusSet_um', 'TargetFillFactor', 'MinDistanceFactor', ...
        'Seed', 'Objective', 'PhaseHeight_nm'});
    resultPath = fullfile(outputDir, 'multi_radius_phase_doe_results.mat');
    bestDesign = struct();
    hazeRiskLimit = 0.12;
    doePreset = "fixture";
    save(resultPath, 'results', 'bestDesign', 'hazeRiskLimit', 'doePreset');
end

function selectContinuousAperture(scriptPath)
    text = strrep(fileread(scriptPath), 'Gratingflag = 1;', 'Gratingflag = 2;');
    file = fopen(scriptPath, 'w', 'n', 'UTF-8');
    assert(file >= 0, 'Cannot write temporary script.');
    closeFile = onCleanup(@() fclose(file));
    fprintf(file, '%s', text);
end
