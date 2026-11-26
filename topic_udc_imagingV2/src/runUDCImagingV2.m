function [optics, images] = runUDCImagingV2(opticsConfig, layerStack, layerImageDir, sceneImage, encoding)
%runUDCImagingV2 用共用远场积分计算透射核，并在线性传感器平面成像。
    arguments
        opticsConfig (1,1) struct
        layerStack (1,1) struct
        layerImageDir (1,1) string
        sceneImage {mustBeNumericOrLogical} = []
        encoding (1,1) string {mustBeMember(encoding,["linear","srgb"])} = "linear"
    end
    moduleDir = fileparts(fileparts(mfilename('fullpath')));
    root = fileparts(moduleDir);
    savedPath = path;
    restorePath = onCleanup(@() path(savedPath));
    addpath(fullfile(root,'topic_coe_reflection','src'));
    [layerGeometry, regionMap, regionPatterns] = displayoptics.readLayerGeometry( ...
        opticsConfig,layerStack,layerImageDir,true);
    validateattributes(opticsConfig.wavelengths, {'numeric'}, {'row','real','finite','positive','nonempty'});
    validateattributes(opticsConfig.focalLength, {'numeric'}, {'scalar','real','finite','positive'});
    validateattributes(opticsConfig.sensorPitch, {'numeric'}, {'scalar','real','finite','positive'});
    validateattributes(opticsConfig.outputSize, {'numeric'}, {'scalar','integer','>=',3,'finite'});
    validateattributes(opticsConfig.receiverBinSamples, {'numeric'}, {'scalar','integer','positive','finite'});
    validateattributes(opticsConfig.pupilRadius, {'numeric'}, {'scalar','real','positive','nonnan'});
    validateattributes(opticsConfig.heightMap, {'numeric'}, {'2d','real','finite'});
    validateattributes(opticsConfig.deltaIndex, {'numeric'}, {'scalar','real','finite'});
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
    sensorAngles = sensorCoordinate/opticsConfig.focalLength;
    sensorEdgeAngle = (sensorHalfWidthPixels+0.5)*opticsConfig.sensorPitch/opticsConfig.focalLength;
    assert(sensorEdgeAngle <= min(opticsConfig.wavelengths)/(2*objectGridStep)*(1+1e-12), ...
        'UDCImagingV2:Aliasing', '传感器窗口超出物面网格的无混叠角域。');

    wavelengthCount = numel(opticsConfig.wavelengths);
    validateattributes(opticsConfig.spectralWeights, {'numeric'}, {'vector','numel',wavelengthCount,'real','finite','nonnegative'});
    validateattributes(opticsConfig.channelResponse, {'numeric'}, {'2d','nrows',wavelengthCount,'real','finite','nonnegative','nonempty'});
    channelSpectralWeights = opticsConfig.channelResponse.*opticsConfig.spectralWeights(:);
    assert(all(sum(channelSpectralWeights,1)>0), 'UDCImagingV2:Spectrum', '各通道均需正的积分光谱权重。');
    channelSpectralWeights = channelSpectralWeights./sum(channelSpectralWeights,1);
    spectralEnergy = zeros(opticsConfig.outputSize,opticsConfig.outputSize,wavelengthCount);
    channelEnergyKernel = zeros(opticsConfig.outputSize,opticsConfig.outputSize,size(channelSpectralWeights,2));
    clearChannelEnergyKernel = zeros(size(channelEnergyKernel));
    [throughput,capturedThroughput,clearCapturedFraction,samplesPerBinAxis] = deal(zeros(wavelengthCount,1));
    receiverMethod = "midpoint";
    if isfield(opticsConfig,'receiverMethod'), receiverMethod = string(opticsConfig.receiverMethod); end
    mustBeMember(receiverMethod,["midpoint","fftBins"]);
    if receiverMethod=="fftBins"
        fftReceiver = corePrepareFftReceiver(opticsConfig,pupil,regionMap,objectGridStep,sensorCoordinate);
        samplesPerBinAxis(:) = NaN; % FFT 分箱使用原生格重叠，无固定子格点数。
    end
    for wavelengthIndex = 1:wavelengthCount
        wavelength = opticsConfig.wavelengths(wavelengthIndex);
        if receiverMethod=="fftBins"
            [spectralEnergy(:,:,wavelengthIndex),clearEnergy,throughput(wavelengthIndex)] = coreIntegrateFftReceiver( ...
                fftReceiver,transmissionAmplitudeByRegion(:,wavelengthIndex),wavelength,sensorCoordinate,opticsConfig.focalLength);
            if wavelengthIndex==1 || mod(wavelengthIndex,10)==0 || wavelengthIndex==wavelengthCount
                fprintf('FFT 接收格积分 %d/%d，%.0f nm\n',wavelengthIndex,wavelengthCount,wavelength*1e9);
            end
        else
            localTransmissionAmplitude = reshape( ...
                transmissionAmplitudeByRegion(double(regionMap(:))+1,wavelengthIndex),objectRows,objectColumns);
            transmittedField = pupil.*localTransmissionAmplitude ...
                .*exp(2i*pi*opticsConfig.deltaIndex*opticsConfig.heightMap/wavelength);
            [spectralEnergy(:,:,wavelengthIndex),~,samplesPerBinAxis(wavelengthIndex)] = displayoptics.integrateFarField( ...
                transmittedField,x,y,sensorAngles,wavelength,opticsConfig.focalLength,opticsConfig.receiverBinSamples,incidentArea);
            clearEnergy = displayoptics.integrateFarField( ...
                pupil,x,y,sensorAngles,wavelength,opticsConfig.focalLength,opticsConfig.receiverBinSamples,incidentArea);
            throughput(wavelengthIndex) = sum(abs(transmittedField).^2,'all')*objectGridStep^2/incidentArea;
        end
        capturedThroughput(wavelengthIndex) = sum(spectralEnergy(:,:,wavelengthIndex),'all');
        clearCapturedFraction(wavelengthIndex) = sum(clearEnergy,'all');
        channelWeightAtWavelength = reshape(channelSpectralWeights(wavelengthIndex,:),1,1,[]);
        channelEnergyKernel = channelEnergyKernel + spectralEnergy(:,:,wavelengthIndex).*channelWeightAtWavelength;
        clearChannelEnergyKernel = clearChannelEnergyKernel + clearEnergy.*channelWeightAtWavelength;
    end
    clear fftReceiver

    optics.config = opticsConfig;
    optics.config.mode = "udc"; % 复用现有线性图像卷积接口。
    optics.config.receiverMethod = receiverMethod;
    optics.layerGeometry = layerGeometry;
    optics.regionMap = regionMap;
    optics.regionPatterns = regionPatterns;
    optics.transmissionLookup = transmissionAmplitudeByRegion;
    optics.pupil = pupil;
    optics.x = x;
    optics.y = y;
    optics.coordinate = sensorCoordinate;
    optics.coordinateUnit = "m (sensor)";
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
    optics.spatialFrequency = (-sensorHalfWidthPixels:sensorHalfWidthPixels) ...
        /(opticsConfig.outputSize*opticsConfig.sensorPitch*1e3); % cycles/mm
    optics.metrics = table(opticsConfig.wavelengths(:)*1e9, ...
        throughput,capturedThroughput,clearCapturedFraction,samplesPerBinAxis, ...
        'VariableNames',{'WavelengthNm','Throughput','CapturedThroughput','ClearCapturedFraction','SamplesPerBinAxis'});
    optics.quality = coreAssessDiffraction(optics);
    images = [];
    if ~isempty(sceneImage)
        images = coreFormUDCScene(sceneImage,optics,encoding);
    end
end
