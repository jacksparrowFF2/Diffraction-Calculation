function transmissionModes = coreCanonicalUDCModeNames(transmissionModes)
%coreCanonicalUDCModeNames 将旧式模式名转换为当前配置名称。
    transmissionModes(transmissionModes == "direct") = "specifiedTransmittance";
end
