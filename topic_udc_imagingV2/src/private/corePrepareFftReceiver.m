function receiver = corePrepareFftReceiver(cfg,pupil,regionMap,dx,coordinate)
%corePrepareFftReceiver 预计算与波长无关的入瞳傅里叶基底。
    assert(isscalar(cfg.heightMap),'UDCImagingV2:FftHeightMap', ...
        'FFT 基底积分需要空间均匀的 heightMap。');
    validateattributes(cfg.receiverFftPaddingFactor,{'numeric'},{'scalar','integer','>=',1,'finite'});
    nfft = 2*ceil(max(size(pupil))*cfg.receiverFftPaddingFactor/2);
    step = 1/(nfft*dx);
    sensorEdge = max(abs(coordinate))+(coordinate(2)-coordinate(1))/2;
    frequencyEdge = sensorEdge/(min(cfg.wavelengths)*cfg.focalLength)+step/2;
    maxBin = max(1,ceil(frequencyEdge/step));
    bins = max(-nfft/2,-maxBin):min(nfft/2-1,maxBin);
    indices = mod(bins,nfft)+1;
    counts = accumarray(double(regionMap(pupil~=0))+1,1,[double(max(regionMap,[],'all'))+1,1]);
    active = find(counts>0);
    [~,order] = sort(counts(active));
    receiver.regions = active(order);
    receiver.clear = fftBand(pupil,nfft,indices);
    receiver.basis = cell(1,numel(active)-1);
    for j = 1:numel(receiver.basis)
        indicator = (regionMap==receiver.regions(j)-1) & pupil~=0;
        receiver.basis{j} = fftBand(indicator,nfft,indices);
    end
    receiver.frequency = bins*step;
    receiver.energyDenominator = nfft^2*nnz(pupil);
    receiver.regionFractions = counts(receiver.regions)/nnz(pupil);
end

function spectrum = fftBand(field,nfft,indices)
%fftBand 分块做两次一维 FFT，每次只保存目标频带。
    blockSize = 128;
    intermediate = complex(zeros(numel(indices),size(field,2),'single'));
    for first = 1:blockSize:size(field,2)
        columns = first:min(first+blockSize-1,size(field,2));
        transformed = fft(single(field(:,columns)),nfft,1);
        intermediate(:,columns) = transformed(indices,:);
    end
    intermediate = intermediate.';
    spectrum = complex(zeros(numel(indices),'single'));
    for first = 1:blockSize:numel(indices)
        columns = first:min(first+blockSize-1,numel(indices));
        transformed = fft(intermediate(:,columns),nfft,1);
        spectrum(:,columns) = transformed(indices,:);
    end
    spectrum = spectrum.';
end
