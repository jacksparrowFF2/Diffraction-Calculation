function modeOutputDir = coreFindUDCModeResults(projectOutputDir,transmissionMode)
%coreFindUDCModeResults 优先读取当前目录名，兼容历史结果目录。
    modeOutputDir = fullfile(projectOutputDir,transmissionMode);
    if transmissionMode == "specifiedTransmittance" && ~isfile(fullfile(modeOutputDir,'udc_results.mat'))
        legacyOutputDir = fullfile(projectOutputDir,'direct');
        if isfile(fullfile(legacyOutputDir,'udc_results.mat'))
            modeOutputDir = legacyOutputDir;
        end
    end
end
