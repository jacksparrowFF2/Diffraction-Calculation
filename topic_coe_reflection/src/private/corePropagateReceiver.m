function [energy, irradianceRatio, samplesPerBin] = corePropagateReceiver(field, x, y, coordinate, lambda, cfg, incidentArea)
%corePropagateReceiver 夫琅禾费接收面强度；不计算观察面二次相位。
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
    % 在指定均匀频率网格上积分；两轴分别卷积，面积元为 dx*dy。
    frequency = fineAngle/lambda;
    spectrum = transformColumns(field.', x, frequency).';
    spectrum = transformColumns(spectrum, y, frequency);
    z = cfg.propagationDistance;
    propagated = spectrum*(exp(2i*pi*z/lambda)/(1i*lambda*z));
    fineIntensity = abs(propagated).^2;
    % MATLAB 行列均按细采样分组，先求每个输出像素内的平均辐照度。
    grouped = reshape(fineIntensity, samplesPerBin, count, samplesPerBin, count);
    irradianceRatio = reshape(sum(sum(grouped, 1), 3), count, count)/samplesPerBin^2;
    energy = irradianceRatio*(z*angularStep)^2/incidentArea;
end

function output = transformColumns(input, coordinate, frequency)
%transformColumns 二次相位卷积和 FFT 加速的一维傅里叶积分，两轴共用。
    count = size(input, 1);
    outputCount = numel(frequency);
    step = coordinate(2) - coordinate(1);
    frequencyStep = frequency(2) - frequency(1);
    index = (0:count-1).';
    outIndex = (0:outputCount-1).';
    difference = [0:outputCount-1, -count:-1].';
    chirped = input.*exp(-1i*pi*(step*frequencyStep*index.^2 ...
        + 2*step*frequency(1)*index));
    kernel = exp(1i*pi*step*frequencyStep*difference.^2);
    convolution = ifft(fft(chirped, count + outputCount, 1).*fft(kernel), [], 1);
    output = convolution(1:outputCount, :).*exp(-1i*pi*( ...
        step*frequencyStep*outIndex.^2 + 2*frequencyStep*coordinate(1)*outIndex ...
        + 2*coordinate(1)*frequency(1)))*step;
end
