function [energy,clearEnergy,throughput] = coreIntegrateFftReceiver(receiver,lookup,lambda,coordinate,focalLength)
%coreIntegrateFftReceiver 在传感器坐标中按格面积重叠积分 FFT 强度。
    native = lambda*focalLength*receiver.frequency;
    nativeStep = native(2)-native(1);
    sensorStep = coordinate(2)-coordinate(1);
    edge = max(abs(coordinate))+(nativeStep+sensorStep)/2;
    active = find(abs(native)<=edge);
    left = max(coordinate(:)-sensorStep/2,native(active)-nativeStep/2);
    right = min(coordinate(:)+sensorStep/2,native(active)+nativeStep/2);
    weights = sparse(max(0,right-left)/nativeStep);
    clear left right

    coefficients = lookup(receiver.regions);
    throughput = sum(receiver.regionFractions.*abs(coefficients).^2);
    reference = single(coefficients(end));
    offsets = single(coefficients(1:end-1)-coefficients(end));
    projected = zeros(numel(coordinate),numel(active));
    projectedClear = projected;
    blockSize = 128;
    for first = 1:blockSize:numel(active)
        block = first:min(first+blockSize-1,numel(active));
        cols = active(block);
        clearSpectrum = receiver.clear(active,cols);
        spectrum = reference*clearSpectrum;
        for j = 1:numel(offsets)
            spectrum = spectrum + offsets(j)*receiver.basis{j}(active,cols);
        end
        projected(:,block) = weights*double(abs(spectrum).^2);
        projectedClear(:,block) = weights*double(abs(clearSpectrum).^2);
    end
    energy = (projected*weights.')/receiver.energyDenominator;
    clearEnergy = (projectedClear*weights.')/receiver.energyDenominator;
end
