function [optics, images] = runUDCImagingV2(cfg, stack, inputDir, scene, encoding)
%runUDCImagingV2 用共用远场积分计算透射核，并在线性传感器平面成像。
    arguments
        cfg (1,1) struct
        stack (1,1) struct
        inputDir (1,1) string
        scene {mustBeNumericOrLogical} = []
        encoding (1,1) string {mustBeMember(encoding,["linear","srgb"])} = "linear"
    end
    moduleDir = fileparts(fileparts(mfilename('fullpath')));
    root = fileparts(moduleDir);
    savedPath = path;
    restorePath = onCleanup(@() path(savedPath));
    addpath(fullfile(root,'topic_coe_reflection','src'), fullfile(root,'topic_udc_imaging','src'));
    [geometry, regionMap, patterns] = displayoptics.readLayerGeometry(cfg, stack, inputDir, true);
    validateattributes(cfg.wavelengths, {'numeric'}, {'row','real','finite','positive','nonempty'});
    validateattributes(cfg.focalLength, {'numeric'}, {'scalar','real','finite','positive'});
    validateattributes(cfg.sensorPitch, {'numeric'}, {'scalar','real','finite','positive'});
    validateattributes(cfg.outputSize, {'numeric'}, {'scalar','integer','>=',3,'finite'});
    validateattributes(cfg.receiverBinSamples, {'numeric'}, {'scalar','integer','positive','finite'});
    validateattributes(cfg.pupilRadius, {'numeric'}, {'scalar','real','positive','nonnan'});
    validateattributes(cfg.heightMap, {'numeric'}, {'2d','real','finite'});
    validateattributes(cfg.deltaIndex, {'numeric'}, {'scalar','real','finite'});
    assert(mod(cfg.outputSize,2)==1, 'UDCImagingV2:OutputSize', '输出采样点数须为奇数。');
    assert(isscalar(cfg.heightMap) || isequal(size(cfg.heightMap),size(regionMap)), ...
        'UDCImagingV2:HeightSize', 'heightMap 须为标量或与计算网格同尺寸。');

    lookup = coreTransmitLayerStack(cfg.wavelengths, stack, patterns, inputDir);
    dx = geometry.gridStep;
    [ny,nx] = size(regionMap);
    x = ((0:nx-1)-floor(nx/2))*dx;
    y = ((0:ny-1)-floor(ny/2))*dx;
    assert(isinf(cfg.pupilRadius) || 2*cfg.pupilRadius <= min(nx,ny)*dx, ...
        'UDCImagingV2:PupilExtent', '有限入瞳必须完整落在计算物面内。');
    pupil = double(y(:).^2+x.^2 <= cfg.pupilRadius^2);
    incidentArea = nnz(pupil)*dx^2;
    assert(incidentArea>0, 'UDCImagingV2:EmptyPupil', '采样后的入瞳为空。');
    half = (cfg.outputSize-1)/2;
    coordinate = (-half:half)*cfg.sensorPitch;
    angles = coordinate/cfg.focalLength;
    edgeAngle = (half+0.5)*cfg.sensorPitch/cfg.focalLength;
    assert(edgeAngle <= min(cfg.wavelengths)/(2*dx)*(1+1e-12), ...
        'UDCImagingV2:Aliasing', '传感器窗口超出物面网格的无混叠角域。');

    nLambda = numel(cfg.wavelengths);
    validateattributes(cfg.spectralWeights, {'numeric'}, {'vector','numel',nLambda,'real','finite','nonnegative'});
    validateattributes(cfg.channelResponse, {'numeric'}, {'2d','nrows',nLambda,'real','finite','nonnegative','nonempty'});
    response = cfg.channelResponse.*cfg.spectralWeights(:);
    assert(all(sum(response,1)>0), 'UDCImagingV2:Spectrum', '各通道均需正的积分光谱权重。');
    response = response./sum(response,1);
    energy = zeros(cfg.outputSize,cfg.outputSize,nLambda);
    kernel = zeros(cfg.outputSize,cfg.outputSize,size(response,2));
    clearKernel = zeros(size(kernel));
    [throughput,captured,clearCaptured,samples] = deal(zeros(nLambda,1));
    method = "midpoint";
    if isfield(cfg,'receiverMethod'), method = string(cfg.receiverMethod); end
    mustBeMember(method,["midpoint","fftBins"]);
    if method=="fftBins"
        fftReceiver = corePrepareFftReceiver(cfg,pupil,regionMap,dx,coordinate);
        samples(:) = NaN; % FFT 分箱使用原生格重叠，无固定子格点数。
    end
    for k = 1:nLambda
        lambda = cfg.wavelengths(k);
        if method=="fftBins"
            [energy(:,:,k),clearEnergy,throughput(k)] = coreIntegrateFftReceiver( ...
                fftReceiver,lookup(:,k),lambda,coordinate,cfg.focalLength);
            if k==1 || mod(k,10)==0 || k==nLambda
                fprintf('FFT 接收格积分 %d/%d，%.0f nm\n',k,nLambda,lambda*1e9);
            end
        else
            t = reshape(lookup(double(regionMap(:))+1,k),ny,nx);
            field = pupil.*t.*exp(2i*pi*cfg.deltaIndex*cfg.heightMap/lambda);
            [energy(:,:,k),~,samples(k)] = displayoptics.integrateFarField( ...
                field,x,y,angles,lambda,cfg.focalLength,cfg.receiverBinSamples,incidentArea);
            clearEnergy = displayoptics.integrateFarField( ...
                pupil,x,y,angles,lambda,cfg.focalLength,cfg.receiverBinSamples,incidentArea);
            throughput(k) = sum(abs(field).^2,'all')*dx^2/incidentArea;
        end
        captured(k) = sum(energy(:,:,k),'all');
        clearCaptured(k) = sum(clearEnergy,'all');
        weight = reshape(response(k,:),1,1,[]);
        kernel = kernel + energy(:,:,k).*weight;
        clearKernel = clearKernel + clearEnergy.*weight;
    end
    clear fftReceiver

    optics.config = cfg;
    optics.config.mode = "udc"; % 复用现有线性图像卷积接口。
    optics.config.receiverMethod = method;
    optics.layerGeometry = geometry;
    optics.regionMap = regionMap;
    optics.regionPatterns = patterns;
    optics.transmissionLookup = lookup;
    optics.pupil = pupil;
    optics.x = x;
    optics.y = y;
    optics.coordinate = coordinate;
    optics.coordinateUnit = "m (sensor)";
    optics.spectralEnergy = energy;
    optics.kernel = kernel;
    optics.clearKernel = clearKernel;
    totals = sum(sum(kernel,1),2);
    optics.channelCapturedThroughput = reshape(totals,1,[]);
    optics.channelThroughput = throughput.'*response;
    totals(totals==0) = 1;
    optics.psf = kernel./totals;
    centeredPsf = ifftshift(ifftshift(optics.psf,1),2);
    otf = fft2(centeredPsf);
    optics.mtf = abs(fftshift(fftshift(otf,1),2));
    optics.spatialFrequency = (-half:half)/(cfg.outputSize*cfg.sensorPitch*1e3); % cycles/mm
    optics.metrics = table(cfg.wavelengths(:)*1e9, throughput, captured, clearCaptured, samples, ...
        'VariableNames',{'WavelengthNm','Throughput','CapturedThroughput','ClearCapturedFraction','SamplesPerBinAxis'});
    optics.quality = coreAssessDiffraction(optics);
    images = [];
    if ~isempty(scene)
        images = formUDCImage(scene,optics,encoding);
    end
end
