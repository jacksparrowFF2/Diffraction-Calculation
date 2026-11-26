function [projectSettings,projectConfigJson,transmissionModes] = coreLoadUDCProjectConfig(projectInputDir,transmissionModes)
%coreLoadUDCProjectConfig 读取配置并按配置顺序选择透射模式。
    projectConfigJson = fileread(fullfile(projectInputDir,'config.json'));
    projectSettings = coreNormalizeUDCModeConfig(jsondecode(projectConfigJson));
    if isempty(transmissionModes)
        transmissionModes = reshape(string(fieldnames(projectSettings.transmissionModes)),1,[]);
    else
        transmissionModes = coreCanonicalUDCModeNames(transmissionModes);
    end
    assert(~isempty(transmissionModes) && all(isfield(projectSettings.transmissionModes,transmissionModes)), ...
        'UDCImagingV2:TransmissionMode','项目须包含所有选定的透射模式。');
end
