function [energy, irradianceRatio, samplesPerBin] = propagateCoeReceiver(field, x, y, coordinate, lambda, cfg, incidentArea)
%propagateCoeReceiver 夫琅禾费接收面强度；不计算观察面二次相位。
% U=exp(i*k*z)/(i*lambda*z)*F，I/Iinc=abs(U)^2。
% 自动细分输出像素，以 midpoint 积分保留窄峰能量；不事后归一化捕获能量。
    angularStep = coordinate(2) - coordinate(1);
    extent = max([numel(x)*(x(2)-x(1)), numel(y)*(y(2)-y(1))]);
    % 每个 lambda/L 至少四个积分采样；奇数细分保持像素中心在采样点上。
    samplesPerBin = max(cfg.receiverBinSamples, ceil(4*extent*angularStep/lambda));
    samplesPerBin = samplesPerBin + mod(samplesPerBin + 1, 2);
    count = numel(coordinate);
    fineCount = count*samplesPerBin;
    fineStep = angularStep/samplesPerBin;
    fineAngle = coordinate(1) - angularStep/2 + ((0:fineCount-1) + 0.5)*fineStep;
    spectrum = fourierIntegralCoe(field, x, y, fineAngle/lambda, fineAngle/lambda);
    z = cfg.propagationDistance;
    propagated = spectrum*(exp(2i*pi*z/lambda)/(1i*lambda*z));
    fineIntensity = abs(propagated).^2;
    % MATLAB 行列均按细采样分组，先求每个输出像素内的平均辐照度。
    grouped = reshape(fineIntensity, samplesPerBin, count, samplesPerBin, count);
    irradianceRatio = reshape(sum(sum(grouped, 1), 3), count, count)/samplesPerBin^2;
    energy = irradianceRatio*(z*angularStep)^2/incidentArea;
end
