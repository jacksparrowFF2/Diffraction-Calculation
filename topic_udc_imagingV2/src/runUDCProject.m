function [optics,images,reportPath] = runUDCProject(projectName,transmissionModes)
%runUDCProject 按项目配置运行选定模式并输出报告。
    arguments
        projectName (1,1) string
        transmissionModes (1,:) string = strings(1,0)
    end
    assert(isvarname(char(projectName)), 'UDCImagingV2:ProjectName','项目名须为有效 MATLAB 标识符。');
    moduleDir = fileparts(fileparts(mfilename('fullpath')));
    repositoryDir = fileparts(moduleDir);
    projectInputDir = fullfile(moduleDir,'Input',projectName);
    [projectSettings,projectConfigJson,transmissionModes] = coreLoadUDCProjectConfig(projectInputDir,transmissionModes);
    if isfield(projectSettings,'reuseProject')
        [optics,images,reportPath] = coreReuseUDCProject(moduleDir,projectName,projectSettings,projectConfigJson,transmissionModes);
        return
    end
    [opticsConfig,layerStack,layerImageDir,sceneImage] = corePrepareUDCProjectInputs(projectSettings,repositoryDir);
    opticsConfig = corePrepareUDCSpectrum(opticsConfig,projectSettings,repositoryDir);

    %% 3. 计算并分别导出各模式的 PSF、MTF 和成像结果
    for transmissionMode = transmissionModes
        layerStack.layers = projectSettings.transmissionModes.(transmissionMode);
        modeOutputDir = fullfile(moduleDir,'Output',projectName,transmissionMode);
        fprintf('UDC V2 %s / %s: %d wavelengths, %d x %d sensor bins\n', ...
            projectName,transmissionMode,numel(opticsConfig.wavelengths),opticsConfig.outputSize,opticsConfig.outputSize);
        [optics,images] = runUDCImagingV2(opticsConfig, layerStack, layerImageDir, sceneImage, "srgb");
        coreSaveUDCResults(modeOutputDir,optics,images,layerStack,projectConfigJson,layerImageDir,transmissionMode);
    end
    reportPath = writeUDCReport(fullfile(moduleDir,'Output',projectName),transmissionModes);
    fprintf('HTML report: %s\n',reportPath);
end
