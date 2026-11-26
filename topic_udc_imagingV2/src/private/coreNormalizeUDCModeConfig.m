function projectSettings = coreNormalizeUDCModeConfig(projectSettings)
%coreNormalizeUDCModeConfig 统一新旧透射模式名称，并保留配置顺序。
    modeSettings = projectSettings.transmissionModes;
    configuredNames = string(fieldnames(modeSettings));
    canonicalNames = coreCanonicalUDCModeNames(configuredNames);
    if ~isequal(configuredNames,canonicalNames)
        assert(numel(unique(canonicalNames)) == numel(canonicalNames), ...
            'UDCImagingV2:TransmissionMode','同一配置不能同时定义 direct 和 specifiedTransmittance。');
        projectSettings.transmissionModes = cell2struct(struct2cell(modeSettings),cellstr(canonicalNames),1);
    end
end
