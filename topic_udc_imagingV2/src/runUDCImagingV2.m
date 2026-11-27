function [optics, images] = runUDCImagingV2(opticsConfig, layerStack, layerImageDir, sceneImage, encoding, sceneModes)
%runUDCImagingV2 用共用远场积分计算透射核，并在线性传感器平面成像。
% 指定 sceneModes 时，返回与模式顺序相同的结果，并共享入瞳 FFT 基底。
    arguments
        opticsConfig (1,1) struct
        layerStack (1,1) struct
        layerImageDir (1,1) string
        sceneImage {mustBeNumericOrLogical} = []
        encoding (1,1) string {mustBeMember(encoding,["linear","srgb"])} = "linear"
        sceneModes (1,:) string = strings(1,0)
    end
    moduleDir = fileparts(fileparts(mfilename('fullpath')));
    root = fileparts(moduleDir);
    savedPath = path;
    restorePath = onCleanup(@() path(savedPath));
    addpath(fullfile(root,'topic_coe_reflection','src'));
    validateattributes(opticsConfig.wavelengths, {'numeric'}, {'row','real','finite','positive','nonempty'});
    validateattributes(opticsConfig.focalLength, {'numeric'}, {'scalar','real','finite','positive'});
    validateattributes(opticsConfig.sensorPitch, {'numeric'}, {'scalar','real','finite','positive'});
    validateattributes(opticsConfig.outputSize, {'numeric'}, {'scalar','integer','>=',3,'finite'});
    validateattributes(opticsConfig.receiverBinSamples, {'numeric'}, {'scalar','integer','positive','finite'});
    validateattributes(opticsConfig.pupilRadius, {'numeric'}, {'scalar','real','positive','nonnan'});
    validateattributes(opticsConfig.heightMap, {'numeric'}, {'2d','real','finite'});
    validateattributes(opticsConfig.deltaIndex, {'numeric'}, {'scalar','real','finite'});
    if isempty(sceneModes)
        sceneModes = "";
    else
        mustBeMember(sceneModes,["sensor","object"]);
    end
    sceneConfigs = cell(size(sceneModes));
    imageDistances = zeros(size(sceneModes));
    for sceneIndex = 1:numel(sceneModes)
        sceneConfig = opticsConfig;
        if sceneModes(sceneIndex)~="", sceneConfig.sceneGeometry.mode = sceneModes(sceneIndex); end
        [sceneConfigs{sceneIndex},imageDistances(sceneIndex)] = resolveScene(sceneConfig);
    end

    [layerGeometry, regionMap, regionPatterns] = displayoptics.readLayerGeometry( ...
        opticsConfig,layerStack,layerImageDir,true);
    assert(mod(opticsConfig.outputSize,2)==1, 'UDCImagingV2:OutputSize', '输出采样点数须为奇数。');
    assert(isscalar(opticsConfig.heightMap) || isequal(size(opticsConfig.heightMap),size(regionMap)), ...
        'UDCImagingV2:HeightSize', 'heightMap 须为标量或与计算网格同尺寸。');

    transmissionAmplitudeByRegion = coreTransmitLayerStack( ...
        opticsConfig.wavelengths,layerStack,regionPatterns,layerImageDir);
    objectGridStep = layerGeometry.gridStep;
    [objectRows,objectColumns] = size(regionMap);
    x = ((0:objectColumns-1)-floor(objectColumns/2))*objectGridStep;
    y = ((0:objectRows-1)-floor(objectRows/2))*objectGridStep;
    assert(isinf(opticsConfig.pupilRadius) || 2*opticsConfig.pupilRadius <= min(objectColumns,objectRows)*objectGridStep, ...
        'UDCImagingV2:PupilExtent', '有限入瞳必须完整落在计算物面内。');
    pupil = double(y(:).^2+x.^2 <= opticsConfig.pupilRadius^2);
    incidentArea = nnz(pupil)*objectGridStep^2;
    assert(incidentArea>0, 'UDCImagingV2:EmptyPupil', '采样后的入瞳为空。');
    sensorHalfWidthPixels = (opticsConfig.outputSize-1)/2;
    sensorCoordinate = (-sensorHalfWidthPixels:sensorHalfWidthPixels)*opticsConfig.sensorPitch;
    sensorEdgeAngle = (sensorHalfWidthPixels+0.5)*opticsConfig.sensorPitch/min(imageDistances);
    assert(sensorEdgeAngle <= min(opticsConfig.wavelengths)/(2*objectGridStep)*(1+1e-12), ...
        'UDCImagingV2:Aliasing', '传感器窗口超出物面网格的无混叠角域。');

    wavelengthCount = numel(opticsConfig.wavelengths);
    validateattributes(opticsConfig.spectralWeights, {'numeric'}, {'vector','numel',wavelengthCount,'real','finite','nonnegative'});
    validateattributes(opticsConfig.channelResponse, {'numeric'}, {'2d','nrows',wavelengthCount,'real','finite','nonnegative'});
    channelSpectralWeights = opticsConfig.channelResponse.*opticsConfig.spectralWeights(:);
    assert(all(sum(channelSpectralWeights,1)>0), 'UDCImagingV2:Spectrum', '各通道均需正的积分光谱权重。');
    channelSpectralWeights = channelSpectralWeights./sum(channelSpectralWeights,1);

    receiverMethod = "midpoint";
    if isfield(opticsConfig,'receiverMethod'), receiverMethod = string(opticsConfig.receiverMethod); end
    mustBeMember(receiverMethod,["midpoint","fftBins"]);
    fftReceiver = [];
    binBounds = zeros(2,numel(imageDistances));
    if receiverMethod=="fftBins"
        [fftReceiver,binBounds] = corePrepareFftReceiver( ...
            opticsConfig,pupil,regionMap,objectGridStep,sensorCoordinate,imageDistances);
    end

    baseOptics.layerGeometry = layerGeometry;
    baseOptics.regionMap = regionMap;
    baseOptics.regionPatterns = regionPatterns;
    baseOptics.transmissionLookup = transmissionAmplitudeByRegion;
    baseOptics.pupil = pupil;
    baseOptics.x = x;
    baseOptics.y = y;
    baseOptics.coordinate = sensorCoordinate;
    baseOptics.coordinateUnit = "m (sensor)";
    baseOptics.spatialFrequency = (-sensorHalfWidthPixels:sensorHalfWidthPixels) ...
        /(opticsConfig.outputSize*opticsConfig.sensorPitch*1e3);

    opticalResults = cell(size(sceneModes));
    for sceneIndex = 1:numel(sceneModes)
        baseOptics.config = sceneConfigs{sceneIndex};
        baseOptics.config.mode = "udc"; % 复用现有线性图像卷积接口。
        baseOptics.config.receiverMethod = receiverMethod;
        if numel(sceneModes)>1
            fprintf('场景 %s：像距 %.6f mm\n',sceneModes(sceneIndex),1e3*imageDistances(sceneIndex));
        end
        opticalResults{sceneIndex} = integrateScene(baseOptics,channelSpectralWeights, ...
            imageDistances(sceneIndex),incidentArea,fftReceiver,binBounds(:,sceneIndex));
    end
    clear fftReceiver
    optics = [opticalResults{:}];
    images = [];
    if ~isempty(sceneImage)
        sceneImages = cell(size(sceneModes));
        for sceneIndex = 1:numel(sceneModes)
            sceneImages{sceneIndex} = coreFormUDCScene(sceneImage,optics(sceneIndex),encoding);
        end
        images = [sceneImages{:}];
    end
end

function optics = integrateScene(optics,channelSpectralWeights,imageDistance,incidentArea,fftReceiver,binBounds)
    cfg = optics.config;
    sensorCoordinate = optics.coordinate;
    sensorAngles = sensorCoordinate/imageDistance;
    wavelengthCount = numel(cfg.wavelengths);
    spectralEnergy = zeros(cfg.outputSize,cfg.outputSize,wavelengthCount);
    channelEnergyKernel = zeros(cfg.outputSize,cfg.outputSize,size(channelSpectralWeights,2));
    clearChannelEnergyKernel = zeros(size(channelEnergyKernel));
    [throughput,capturedThroughput,clearCapturedFraction,samplesPerBinAxis] = deal(zeros(wavelengthCount,1));
    if cfg.receiverMethod=="fftBins"
        samplesPerBinAxis(:) = NaN; % FFT 分箱无固定子格点数。
    end
    for wavelengthIndex = 1:wavelengthCount
        wavelength = cfg.wavelengths(wavelengthIndex);
        if cfg.receiverMethod=="fftBins"
            [spectralEnergy(:,:,wavelengthIndex),clearEnergy,throughput(wavelengthIndex)] = coreIntegrateFftReceiver( ...
                fftReceiver,optics.transmissionLookup(:,wavelengthIndex),wavelength,sensorCoordinate,imageDistance,binBounds);
            if wavelengthIndex==1 || mod(wavelengthIndex,10)==0 || wavelengthIndex==wavelengthCount
                fprintf('FFT 接收格积分 %d/%d，%.0f nm\n',wavelengthIndex,wavelengthCount,wavelength*1e9);
            end
        else
            localTransmissionAmplitude = reshape( ...
                optics.transmissionLookup(double(optics.regionMap(:))+1,wavelengthIndex),size(optics.regionMap));
            transmittedField = optics.pupil.*localTransmissionAmplitude ...
                .*exp(2i*pi*cfg.deltaIndex*cfg.heightMap/wavelength);
            [spectralEnergy(:,:,wavelengthIndex),~,samplesPerBinAxis(wavelengthIndex)] = displayoptics.integrateFarField( ...
                transmittedField,optics.x,optics.y,sensorAngles,wavelength,imageDistance,cfg.receiverBinSamples,incidentArea);
            clearEnergy = displayoptics.integrateFarField( ...
                optics.pupil,optics.x,optics.y,sensorAngles,wavelength,imageDistance,cfg.receiverBinSamples,incidentArea);
            throughput(wavelengthIndex) = sum(abs(transmittedField).^2,'all')*optics.layerGeometry.gridStep^2/incidentArea;
        end
        capturedThroughput(wavelengthIndex) = sum(spectralEnergy(:,:,wavelengthIndex),'all');
        clearCapturedFraction(wavelengthIndex) = sum(clearEnergy,'all');
        channelWeightAtWavelength = reshape(channelSpectralWeights(wavelengthIndex,:),1,1,[]);
        channelEnergyKernel = channelEnergyKernel + spectralEnergy(:,:,wavelengthIndex).*channelWeightAtWavelength;
        clearChannelEnergyKernel = clearChannelEnergyKernel + clearEnergy.*channelWeightAtWavelength;
    end

    optics.spectralEnergy = spectralEnergy;
    optics.kernel = channelEnergyKernel;
    optics.clearKernel = clearChannelEnergyKernel;
    psfNormalization = sum(sum(channelEnergyKernel,1),2);
    optics.channelCapturedThroughput = reshape(psfNormalization,1,[]);
    optics.channelThroughput = throughput.'*channelSpectralWeights;
    psfNormalization(psfNormalization==0) = 1;
    optics.psf = channelEnergyKernel./psfNormalization;
    centeredPsf = ifftshift(ifftshift(optics.psf,1),2);
    otf = fft2(centeredPsf);
    optics.mtf = abs(fftshift(fftshift(otf,1),2));
    optics.metrics = table(cfg.wavelengths(:)*1e9, ...
        throughput,capturedThroughput,clearCapturedFraction,samplesPerBinAxis, ...
        'VariableNames',{'WavelengthNm','Throughput','CapturedThroughput','ClearCapturedFraction','SamplesPerBinAxis'});
    optics.quality = coreAssessDiffraction(optics);
end

function [opticsConfig,imageDistance] = resolveScene(opticsConfig)
    imageDistance = opticsConfig.focalLength;
    if isfield(opticsConfig,'sceneProjection'), opticsConfig = rmfield(opticsConfig,'sceneProjection'); end
    if isfield(opticsConfig,'sceneGeometry')
        geometry = opticsConfig.sceneGeometry;
        validateattributes(geometry,{'struct'},{'scalar'});
        mode = string(geometry.mode);
        validateattributes(mode,{'string'},{'scalar'});
        mustBeMember(mode,["sensor","object"]);
        if mode=="object"
            validateattributes(geometry.objectDistanceMm,{'numeric'}, ...
                {'scalar','real','finite','>',opticsConfig.focalLength*1e3});
            validateattributes(geometry.objectSizeMm,{'numeric'}, ...
                {'vector','numel',2,'real','finite','positive'});
            objectDistance = geometry.objectDistanceMm*1e-3;
            magnification = opticsConfig.focalLength/(objectDistance-opticsConfig.focalLength);
            % 对焦时，薄透镜的二次相位抵消物面及像面的轴上曲率。
            imageDistance = opticsConfig.focalLength*(1+magnification);
            imageSizeMm = reshape(geometry.objectSizeMm,1,2)*magnification;
            opticsConfig.sceneProjection = struct('imageDistanceMm',1e3*imageDistance, ...
                'magnification',magnification,'imageSizeMm',imageSizeMm, ...
                'imageSizePixels',max(1,round(imageSizeMm*1e-3/opticsConfig.sensorPitch)));
            validateattributes(opticsConfig.sceneProjection.imageSizePixels,{'numeric'},{'finite'});
            if isfield(geometry,'sensorFrameMm')
                validateattributes(geometry.sensorFrameMm,{'numeric'}, ...
                    {'vector','numel',2,'real','finite','positive'});
                framePixels = max(1,round(reshape(geometry.sensorFrameMm,1,2)*1e-3/opticsConfig.sensorPitch));
                validateattributes(framePixels,{'numeric'},{'finite'});
                previewWidth = 1200;
                if isfield(geometry,'previewWidthPixels'), previewWidth = geometry.previewWidthPixels; end
                validateattributes(previewWidth,{'numeric'},{'scalar','real','finite','integer','positive'});
                previewWidth = min(previewWidth,framePixels(1));
                previewHeight = max(1,round((previewWidth/framePixels(1))*framePixels(2)));
                background = 0.18;
                if isfield(geometry,'backgroundLinear'), background = geometry.backgroundLinear; end
                validateattributes(background,{'numeric'},{'scalar','real','finite','>=',0,'<=',1});
                opticsConfig.sceneProjection.frame = struct('sizePixels',framePixels, ...
                    'sizeMm',framePixels*(opticsConfig.sensorPitch*1e3), ...
                    'previewSizePixels',[previewWidth,previewHeight],'backgroundLinear',background);
            end
        end
    end
end
