function [optics,images,reportPath] = runUDCProject(projectName,transmissionModes)
%runUDCProject 按项目配置运行选定模式并输出报告。
    arguments
        projectName (1,1) string
        transmissionModes (1,:) string = strings(1,0)
    end
    assert(isvarname(char(projectName)), 'UDCImagingV2:ProjectName','项目名须为有效 MATLAB 标识符。');
    moduleDir = fileparts(fileparts(mfilename('fullpath')));
    repositoryDir = fileparts(moduleDir);
    inputDir = fullfile(moduleDir,'Input',projectName);
    projectConfig = fileread(fullfile(inputDir,'config.json'));
    settings = jsondecode(projectConfig);
    if isempty(transmissionModes)
        transmissionModes = reshape(string(fieldnames(settings.transmissionModes)),1,[]);
    end
    cfg = settings.cfg;
    stack = settings.stack;
    cfg.periods = reshape(cfg.periods,1,[]);
    stack.layerNames = reshape(string(stack.layerNames),1,[]);
    stack.layerLevels = reshape(stack.layerLevels,1,[]);
    imageDir = fullfile(repositoryDir,string(settings.layerImageDir));
    scene = imread(fullfile(repositoryDir,string(settings.sceneFile)));

    %% 2. D65 梯形积分及示例相机响应
    band = settings.wavelengthRangeNm;
    wavelengthNm = band(1):settings.wavelengthStepNm:band(2);
    assert(numel(wavelengthNm)>=2 && wavelengthNm(end)==band(2), ...
        'UDCImagingV2:WavelengthGrid', '波长步长须整除计算波段。');
    cfg.wavelengths = wavelengthNm*1e-9;
    d65 = readmatrix(fullfile(repositoryDir,'topic_coe_reflection','Input','Colorimetry','CIE_std_illum_D65.csv'));
    cfg.spectralWeights = interp1(d65(:,1),d65(:,2),wavelengthNm,'linear')*settings.wavelengthStepNm;
    cfg.spectralWeights([1,end]) = cfg.spectralWeights([1,end])/2;
    centers = reshape(settings.channelCentersNm,1,[]);
    sigma = reshape(settings.channelSigmaNm,1,[]);
    cfg.channelResponse = exp(-0.5*((wavelengthNm(:)-centers)./sigma).^2);

    %% 3. 计算并分别导出各模式的 PSF、MTF 和成像结果
    for mode = transmissionModes
        stack.layers = settings.transmissionModes.(mode);
        outputDir = fullfile(moduleDir,'Output',projectName,mode);
        fprintf('UDC V2 %s / %s: %d wavelengths, %d x %d sensor bins\n', ...
            projectName,mode,numel(wavelengthNm),cfg.outputSize,cfg.outputSize);
        [optics,images] = runUDCImagingV2(cfg, stack, imageDir, scene, "srgb");
        saveResults(outputDir,optics,images,stack,projectConfig,imageDir,mode);
    end
    reportPath = writeUDCReport(fullfile(moduleDir,'Output',projectName),transmissionModes);
    fprintf('HTML report: %s\n',reportPath);
end

function saveResults(outputDir,optics,images,stack,projectConfig,imageDir,mode)
    cfg = optics.config;
    if ~isfolder(outputDir), mkdir(outputDir); end
    writetable(optics.metrics,fullfile(outputDir,'udc_transmission_spectrum.csv'));
    channelMetrics = withRgbLabels(optics.quality.channelMetrics);
    writetable(channelMetrics,fullfile(outputDir,'udc_channel_metrics.csv'));
    writetable(withRgbLabels(optics.quality.tail),fullfile(outputDir,'udc_psf_tail.csv'));
    writetable(withRgbLabels(optics.quality.mtfProfiles),fullfile(outputDir,'udc_mtf_directional.csv'));
    writetable(withRgbLabels(optics.quality.mtfSummary),fullfile(outputDir,'udc_mtf_directional_summary.csv'));
    center = (cfg.outputSize+1)/2;
    mtf = table(optics.spatialFrequency(center:end).',squeeze(optics.mtf(center,center:end,:)), ...
        'VariableNames',{'CyclesPerMm','HorizontalMTF_RGB'});
    writetable(mtf,fullfile(outputDir,'udc_mtf.csv'));
    save(fullfile(outputDir,'udc_results.mat'),'optics','images','cfg','stack','projectConfig','imageDir','mode','-v7.3');
    imwrite(encodeSrgb(images.clear),fullfile(outputDir,'udc_clear.png'));
    imwrite(encodeSrgb(images.fixedExposure),fullfile(outputDir,'udc_fixed_exposure.png'));
    imwrite(encodeSrgb(images.brightnessMatched),fullfile(outputDir,'udc_brightness_matched.png'));
    plotDiffractionDiagnostics(optics,outputDir);
    plotImagingComparison(optics,outputDir);
    disp(channelMetrics);
    fprintf('Results: %s\n',outputDir);
end

function metrics = withRgbLabels(metrics)
    rgb = ["R";"G";"B"];
    metrics.Channel = rgb(metrics.Channel);
end

function encoded = encodeSrgb(linear)
    linear = min(max(linear,0),1);
    encoded = 12.92*linear;
    high = linear > 0.0031308;
    encoded(high) = 1.055*linear(high).^(1/2.4)-0.055;
end
