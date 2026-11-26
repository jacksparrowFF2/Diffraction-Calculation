function coreSaveUDCResults(outputDir,optics,images,stack,projectConfig,imageDir,mode)
%coreSaveUDCResults 保存线性仿真数据、指标表、场景 PNG 和诊断图。
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
