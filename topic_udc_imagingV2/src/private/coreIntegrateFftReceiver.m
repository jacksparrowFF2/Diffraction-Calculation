function [energy,clearEnergy,throughput] = coreIntegrateFftReceiver(receiver,transmissionAmplitudeByRegion,wavelength,sensorCoordinate,focalLength)
%coreIntegrateFftReceiver 在传感器坐标中按格面积重叠积分 FFT 强度。
    fftBinCoordinates = wavelength*focalLength*receiver.frequency;
    fftBinWidth = fftBinCoordinates(2)-fftBinCoordinates(1);
    sensorBinWidth = sensorCoordinate(2)-sensorCoordinate(1);
    overlapLimit = max(abs(sensorCoordinate))+(fftBinWidth+sensorBinWidth)/2;
    overlappingBins = find(abs(fftBinCoordinates)<=overlapLimit);
    overlapStart = max(sensorCoordinate(:)-sensorBinWidth/2,fftBinCoordinates(overlappingBins)-fftBinWidth/2);
    overlapEnd = min(sensorCoordinate(:)+sensorBinWidth/2,fftBinCoordinates(overlappingBins)+fftBinWidth/2);
    overlapWeights = sparse(max(0,overlapEnd-overlapStart)/fftBinWidth);
    clear overlapStart overlapEnd

    regionAmplitudes = transmissionAmplitudeByRegion(receiver.regions);
    throughput = sum(receiver.regionFractions.*abs(regionAmplitudes).^2);
    referenceAmplitude = single(regionAmplitudes(end));
    amplitudeDifferences = single(regionAmplitudes(1:end-1)-regionAmplitudes(end));
    energyIntegratedAlongRows = zeros(numel(sensorCoordinate),numel(overlappingBins));
    clearEnergyIntegratedAlongRows = energyIntegratedAlongRows;
    blockSize = 128;
    for first = 1:blockSize:numel(overlappingBins)
        block = first:min(first+blockSize-1,numel(overlappingBins));
        fftColumns = overlappingBins(block);
        clearSpectrum = receiver.clear(overlappingBins,fftColumns);
        spectrum = referenceAmplitude*clearSpectrum;
        for j = 1:numel(amplitudeDifferences)
            spectrum = spectrum + amplitudeDifferences(j)*receiver.basis{j}(overlappingBins,fftColumns);
        end
        energyIntegratedAlongRows(:,block) = overlapWeights*double(abs(spectrum).^2);
        clearEnergyIntegratedAlongRows(:,block) = overlapWeights*double(abs(clearSpectrum).^2);
    end
    energy = (energyIntegratedAlongRows*overlapWeights.')/receiver.energyDenominator;
    clearEnergy = (clearEnergyIntegratedAlongRows*overlapWeights.')/receiver.energyDenominator;
end
