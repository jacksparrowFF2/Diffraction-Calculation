%% UDC 透射衍射与成像：输入 → 光谱 → 衍射 → 成像 → 导出

%% 1. 选择项目和透射模式
projectName = "case9";
transmissionModes = strings(1,0); % 留空时运行配置中的全部模式。
% "specifiedTransmittance" 给定功率透过率，"nk" 按薄膜 n、k 与厚度计算。
assert(isvarname(char(projectName)), ...
    'UDCImagingV2:ProjectName','项目名须为有效 MATLAB 标识符。');
moduleDir = fileparts(fileparts(mfilename('fullpath')));
repositoryDir = fileparts(moduleDir);
projectInputDir = fullfile(moduleDir,'Input',projectName);
projectOutputDir = fullfile(moduleDir,'Output',projectName);
[projectSettings,projectConfigJson,transmissionModes] = coreLoadUDCProjectConfig(projectInputDir,transmissionModes);
if isfield(projectSettings,'sceneModes')
    [optics,images,reportPath] = runUDCProject(projectName,transmissionModes);
    return
end

%% 2. 读取模型参数、膜堆配置和场景
if isfield(projectSettings,'reuseProject')
    % 复用已计算的基准项目，按配置施加中性衰减并导出结果。
    [optics,images,reportPath] = coreReuseUDCProject( ...
        moduleDir,projectName,projectSettings,projectConfigJson,transmissionModes);
    return
end
[opticsConfig,layerStack,layerImageDir,sceneImage] = corePrepareUDCProjectInputs(projectSettings,repositoryDir);

%% 3. 建立 D65 光谱积分权重和相机 RGB 响应
opticsConfig = corePrepareUDCSpectrum(opticsConfig,projectSettings,repositoryDir);

%% 4. 逐模式计算、成像并保存
for transmissionMode = transmissionModes
    layerStack.layers = projectSettings.transmissionModes.(transmissionMode);
    modeOutputDir = fullfile(projectOutputDir,transmissionMode);
    fprintf('UDC V2 %s / %s: %d wavelengths, %d x %d sensor bins\n', ...
        projectName,transmissionMode,numel(opticsConfig.wavelengths),opticsConfig.outputSize,opticsConfig.outputSize);

    % 4.1 读取膜层图案，求复振幅透射系数，积分传感器接收格能量并计算 PSF、MTF。
    optics = runUDCImagingV2(opticsConfig,layerStack,layerImageDir);

    % 4.2 在线性传感器平面卷积，得到无遮挡、同曝光和亮度匹配图像。
    images = coreFormUDCScene(sceneImage,optics,"srgb");

    % 4.3 导出 MAT、指标 CSV、场景 PNG 和诊断图。
    coreSaveUDCResults(modeOutputDir,optics,images,layerStack,projectConfigJson,layerImageDir,transmissionMode);
end

%% 5. 生成项目 HTML 报告
reportPath = writeUDCReport(projectOutputDir,transmissionModes);
fprintf('HTML report: %s\n',reportPath);
% optics 和 images 保留最后一个模式的结果；各模式的数据分别保存在对应目录。
