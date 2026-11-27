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
    sceneModes = strings(1,0);
    if isfield(projectSettings,'sceneModes')
        validateattributes(projectSettings.sceneModes,{'cell'},{'vector','nonempty'});
        sceneModes = reshape(string(projectSettings.sceneModes),1,[]);
        mustBeMember(sceneModes,["sensor","object"]);
        assert(numel(unique(sceneModes))==numel(sceneModes),'UDCImagingV2:SceneModes','sceneModes 不能重复。');
        assert(~isfield(projectSettings,'reuseProject'),'UDCImagingV2:SceneModes', ...
            '多场景入口需要计算各模式的场景与对应像距。');
        assert(~any(sceneModes=="object") || isfield(projectSettings.cfg,'sceneGeometry'), ...
            'UDCImagingV2:SceneModes','object 模式须配置 cfg.sceneGeometry。');
    end
    if isfield(projectSettings,'reuseProject')
        [optics,images,reportPath] = coreReuseUDCProject(moduleDir,projectName,projectSettings,projectConfigJson,transmissionModes);
        return
    end
    [opticsConfig,layerStack,layerImageDir,sceneImage] = corePrepareUDCProjectInputs(projectSettings,repositoryDir);
    opticsConfig = corePrepareUDCSpectrum(opticsConfig,projectSettings,repositoryDir);

    %% 3. 计算并分别导出各模式的 PSF、MTF 和成像结果
    projectOutputDir = fullfile(moduleDir,'Output',projectName);
    activeScenes = sceneModes;
    if isempty(activeScenes), activeScenes = ""; end
    runLabel = projectName;
    if ~isempty(sceneModes), runLabel = projectName+" / "+strjoin(sceneModes," + "); end
    for transmissionMode = transmissionModes
        clear opticalResults sceneImages optics images
        layerStack.layers = projectSettings.transmissionModes.(transmissionMode);
        fprintf('UDC V2 %s / %s: %d wavelengths, %d x %d sensor bins\n', ...
            runLabel,transmissionMode,numel(opticsConfig.wavelengths),opticsConfig.outputSize,opticsConfig.outputSize);
        [opticalResults,sceneImages] = runUDCImagingV2(opticsConfig,layerStack,layerImageDir,sceneImage,"srgb",sceneModes);
        for sceneIndex = 1:numel(activeScenes)
            sceneOutputDir = projectOutputDir;
            if activeScenes(sceneIndex)=="object"
                sceneOutputDir = fullfile(projectOutputDir,"object");
            end
            modeOutputDir = fullfile(sceneOutputDir,transmissionMode);
            optics = opticalResults(sceneIndex);
            images = sceneImages(sceneIndex);
            coreSaveUDCResults(modeOutputDir,optics,images,layerStack,projectConfigJson,layerImageDir,transmissionMode);
        end
    end
    reportPath = writeUDCReport(projectOutputDir,transmissionModes,sceneModes);
    fprintf('HTML report: %s\n',reportPath);
end
