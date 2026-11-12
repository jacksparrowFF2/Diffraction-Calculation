classdef ModuleIsolationTest < matlab.unittest.TestCase
    %ModuleIsolationTest Run copied feature folders without sibling dependencies.
    methods (Test)
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

function folder = portableModule(testCase, moduleName)
    root = fileparts(fileparts(mfilename('fullpath')));
    temporary = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture());
    folder = fullfile(temporary.Folder, moduleName);
    mkdir(folder);
    % Never copy historical outputs: output checks must observe newly generated files.
    for part = ["Input", "src", "docs", "dist"]
        copyfile(fullfile(root, moduleName, part), fullfile(folder, part));
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
    cfg.angleHalfRangeDeg = 6;
end

function prepareSmallScript(testCase, scriptPath)
    % Only the temporary copy is changed; production simulation defaults stay intact.
    text = fileread(scriptPath);
    text = strrep(text, 'clear, clc, close all', 'clear');
    text = strrep(text, 'nn = 1000;', 'nn = 24;');
    text = strrep(text, 'nn = 700;', 'nn = 24;');
    text = strrep(text, 'useBestDoeDesign = true;', 'useBestDoeDesign = false;');
    text = strrep(text, 'phaseFillFactor = 0.40;', 'phaseFillFactor = 0.01;');
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
