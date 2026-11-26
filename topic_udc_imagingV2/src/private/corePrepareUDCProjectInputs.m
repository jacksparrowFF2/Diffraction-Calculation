function [opticsConfig,layerStack,layerImageDir,sceneImage] = corePrepareUDCProjectInputs(projectSettings,repositoryDir)
%corePrepareUDCProjectInputs 读取物面、膜堆和场景输入。
    assert(~isfield(projectSettings,'neutralTransmittance'),'UDCImagingV2:NeutralTransmittance', ...
        '中性衰减须指定已有的未衰减基准项目 reuseProject。');
    opticsConfig = projectSettings.cfg;
    layerStack = projectSettings.stack;
    opticsConfig.periods = reshape(opticsConfig.periods,1,[]);
    layerStack.layerNames = reshape(string(layerStack.layerNames),1,[]);
    layerStack.layerLevels = reshape(layerStack.layerLevels,1,[]);
    layerImageDir = fullfile(repositoryDir,string(projectSettings.layerImageDir));
    sceneImage = imread(fullfile(repositoryDir,string(projectSettings.sceneFile)));
end
