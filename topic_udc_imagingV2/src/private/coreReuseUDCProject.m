function [optics,images,reportPath] = coreReuseUDCProject(moduleDir,projectName,projectSettings,projectConfigJson,transmissionModes)
%coreReuseUDCProject 复用基准结果并按需施加均匀中性衰减。
    sourceProjectName = string(projectSettings.reuseProject);
    assert(isscalar(sourceProjectName) && isvarname(char(sourceProjectName)) && sourceProjectName~=projectName, ...
        'UDCImagingV2:ReuseProject', ...
        'reuseProject 须为另一有效项目名。');
    physicalProjectSettings = rmfield(projectSettings,'reuseProject');
    neutralTransmittance = 1;
    if isfield(projectSettings,'neutralTransmittance')
        neutralTransmittance = projectSettings.neutralTransmittance;
        validateattributes(neutralTransmittance,{'numeric'},{'scalar','real','finite','positive','<=',1});
        physicalProjectSettings = rmfield(physicalProjectSettings,'neutralTransmittance');
    end
    outputRoot = fullfile(moduleDir,'Output');
    for transmissionMode = transmissionModes
        sourceModeDir = coreFindUDCModeResults(fullfile(outputRoot,sourceProjectName),transmissionMode);
        sourceResultsFile = fullfile(sourceModeDir,'udc_results.mat');
        assert(isfile(sourceResultsFile),'UDCImagingV2:ReuseProject','找不到基准结果：%s',sourceResultsFile);
        savedResults = load(sourceResultsFile,'optics','images','stack','projectConfig','imageDir');
        sourceProjectSettings = coreNormalizeUDCModeConfig(jsondecode(savedResults.projectConfig));
        assert(~isfield(sourceProjectSettings,'neutralTransmittance') || sourceProjectSettings.neutralTransmittance == 1, ...
            'UDCImagingV2:ReuseSettings','基准结果须为未衰减结果：%s',sourceProjectName);
        if isfield(sourceProjectSettings,'reuseProject')
            sourceProjectSettings = rmfield(sourceProjectSettings,'reuseProject');
        end
        if isfield(sourceProjectSettings,'neutralTransmittance')
            sourceProjectSettings = rmfield(sourceProjectSettings,'neutralTransmittance');
        end
        assert(isequaln(physicalProjectSettings,sourceProjectSettings),'UDCImagingV2:ReuseSettings', ...
            '项目与基准结果的原始配置不一致：%s',sourceProjectName);
        optics = savedResults.optics;
        images = savedResults.images;
        optics.config.reusedFromProject = sourceProjectName;
        optics.config.neutralTransmittance = neutralTransmittance;
        if neutralTransmittance ~= 1
            % 均匀中性衰减只改变场振幅和绝对强度，归一化及无遮挡参考保持不变。
            optics.transmissionLookup = optics.transmissionLookup * sqrt(neutralTransmittance);
            optics.spectralEnergy = optics.spectralEnergy .* neutralTransmittance;
            optics.kernel = optics.kernel .* neutralTransmittance;
            optics.channelThroughput = optics.channelThroughput * neutralTransmittance;
            optics.channelCapturedThroughput = optics.channelCapturedThroughput * neutralTransmittance;
            optics.metrics.Throughput = optics.metrics.Throughput * neutralTransmittance;
            optics.metrics.CapturedThroughput = optics.metrics.CapturedThroughput * neutralTransmittance;
            optics.quality = coreAssessDiffraction(optics);
            images.fixedExposure = images.fixedExposure .* neutralTransmittance;
        end
        fprintf('UDC V2 %s / %s: %s, T = %.3g\n',projectName,transmissionMode,sourceProjectName,neutralTransmittance);
        modeOutputDir = fullfile(outputRoot,projectName,transmissionMode);
        coreSaveUDCResults(modeOutputDir,optics,images,savedResults.stack, ...
            projectConfigJson,savedResults.imageDir,transmissionMode);
        clear savedResults
    end
    reportPath = writeUDCReport(fullfile(outputRoot,projectName),transmissionModes);
    fprintf('HTML report: %s\n',reportPath);
end
