function receiver = corePrepareFftReceiver(opticsConfig,pupil,regionMap,objectGridStep,sensorCoordinate)
%corePrepareFftReceiver 预计算与波长无关的入瞳傅里叶基底。
    assert(isscalar(opticsConfig.heightMap),'UDCImagingV2:FftHeightMap', ...
        'FFT 基底积分需要空间均匀的 heightMap。');
    validateattributes(opticsConfig.receiverFftPaddingFactor,{'numeric'},{'scalar','integer','>=',1,'finite'});
    fftSize = 2*ceil(max(size(pupil))*opticsConfig.receiverFftPaddingFactor/2);
    frequencyStep = 1/(fftSize*objectGridStep);
    sensorEdgeCoordinate = max(abs(sensorCoordinate))+(sensorCoordinate(2)-sensorCoordinate(1))/2;
    requiredFrequencyEdge = sensorEdgeCoordinate/(min(opticsConfig.wavelengths)*opticsConfig.focalLength)+frequencyStep/2;
    maxFrequencyBin = max(1,ceil(requiredFrequencyEdge/frequencyStep));
    frequencyBins = max(-fftSize/2,-maxFrequencyBin):min(fftSize/2-1,maxFrequencyBin);
    fftIndices = mod(frequencyBins,fftSize)+1;
    regionPixelCounts = accumarray(double(regionMap(pupil~=0))+1,1,[double(max(regionMap,[],'all'))+1,1]);
    illuminatedRegions = find(regionPixelCounts>0);
    [~,regionOrder] = sort(regionPixelCounts(illuminatedRegions));
    receiver.regions = illuminatedRegions(regionOrder);
    receiver.clear = fftBand(pupil,fftSize,fftIndices);
    receiver.basis = cell(1,numel(illuminatedRegions)-1);
    for j = 1:numel(receiver.basis)
        regionPupilMask = (regionMap==receiver.regions(j)-1) & pupil~=0;
        receiver.basis{j} = fftBand(regionPupilMask,fftSize,fftIndices);
    end
    receiver.frequency = frequencyBins*frequencyStep;
    receiver.energyDenominator = fftSize^2*nnz(pupil);
    receiver.regionFractions = regionPixelCounts(receiver.regions)/nnz(pupil);
end

function spectrum = fftBand(field,fftSize,fftIndices)
%fftBand 分块做两次一维 FFT，每次只保存目标频带。
    blockSize = 128;
    intermediate = complex(zeros(numel(fftIndices),size(field,2),'single'));
    for first = 1:blockSize:size(field,2)
        columns = first:min(first+blockSize-1,size(field,2));
        transformed = fft(single(field(:,columns)),fftSize,1);
        intermediate(:,columns) = transformed(fftIndices,:);
    end
    intermediate = intermediate.';
    spectrum = complex(zeros(numel(fftIndices),'single'));
    for first = 1:blockSize:numel(fftIndices)
        columns = first:min(first+blockSize-1,numel(fftIndices));
        transformed = fft(intermediate(:,columns),fftSize,1);
        spectrum(:,columns) = transformed(fftIndices,:);
    end
    spectrum = spectrum.';
end
