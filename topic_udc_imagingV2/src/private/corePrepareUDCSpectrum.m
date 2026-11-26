function opticsConfig = corePrepareUDCSpectrum(opticsConfig,projectSettings,repositoryDir)
%corePrepareUDCSpectrum 构建 D65 梯形积分权重及示例相机响应。
    wavelengthRangeNm = projectSettings.wavelengthRangeNm;
    wavelengthNm = wavelengthRangeNm(1):projectSettings.wavelengthStepNm:wavelengthRangeNm(2);
    assert(numel(wavelengthNm)>=2 && wavelengthNm(end)==wavelengthRangeNm(2), ...
        'UDCImagingV2:WavelengthGrid', '波长步长须整除计算波段。');
    opticsConfig.wavelengths = wavelengthNm*1e-9;
    d65Spectrum = readmatrix(fullfile(repositoryDir,'topic_coe_reflection','Input','Colorimetry','CIE_std_illum_D65.csv'));
    opticsConfig.spectralWeights = interp1(d65Spectrum(:,1),d65Spectrum(:,2),wavelengthNm,'linear')*projectSettings.wavelengthStepNm;
    opticsConfig.spectralWeights([1,end]) = opticsConfig.spectralWeights([1,end])/2;
    channelCentersNm = reshape(projectSettings.channelCentersNm,1,[]);
    channelSigmaNm = reshape(projectSettings.channelSigmaNm,1,[]);
    opticsConfig.channelResponse = exp(-0.5*((wavelengthNm(:)-channelCentersNm)./channelSigmaNm).^2);
end
