function [energy, irradianceRatio, samplesPerBin] = corePropagateReceiver(field, x, y, coordinate, lambda, cfg, incidentArea)
%corePropagateReceiver 按 COE 的距离和像素积分设置计算远场接收。
    [energy, irradianceRatio, samplesPerBin] = displayoptics.integrateFarField( ...
        field, x, y, coordinate, lambda, cfg.propagationDistance, cfg.receiverBinSamples, incidentArea);
end
