function spectrum = fourierIntegralCoe(field, x, y, fx, fy)
%fourierIntegralCoe 指定均匀频率网格的傅里叶积分，采用二次相位卷积和 FFT 加速。
% 与 topic_aperture_fraunhofer 的 fourierTransform2D_SPH 相同求和方法；
% 分离执行两轴卷积以降低内存。物面不重复周期终点，面积元为 dx*dy。
    spectrum = transformColumns(field.', x, fx).';
    spectrum = transformColumns(spectrum, y, fy);
end

function output = transformColumns(input, coordinate, frequency)
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
