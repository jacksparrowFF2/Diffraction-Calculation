% Compare case01 at fixed pupil and then at a larger pupil.
benchmarkDir = fileparts(mfilename('fullpath'));
moduleDir = fileparts(fileparts(benchmarkDir));
root = fileparts(moduleDir);
out = fullfile(moduleDir,'Output','case01','array_pupil_comparison');
if ~isfolder(out), mkdir(out); end
addpath(fullfile(moduleDir,'src'));
settings = jsondecode(fileread(fullfile(moduleDir,'Input','case01','config.json')));
cfg = settings.cfg;
stack = settings.stack;
stack.layerNames = reshape(string(stack.layerNames),1,[]);
stack.layerLevels = reshape(stack.layerLevels,1,[]);
stack.layers = settings.transmissionModes.nk;
imageDir = string(fullfile(root,settings.layerImageDir));
scene = imread(fullfile(root,settings.sceneFile));
band = settings.wavelengthRangeNm;
wavelengthNm = band(1):settings.wavelengthStepNm:band(2);
cfg.wavelengths = wavelengthNm*1e-9;
d65 = readmatrix(fullfile(root,'topic_coe_reflection','Input','Colorimetry','CIE_std_illum_D65.csv'));
cfg.spectralWeights = interp1(d65(:,1),d65(:,2),wavelengthNm,'linear')*settings.wavelengthStepNm;
cfg.spectralWeights([1,end]) = cfg.spectralWeights([1,end])/2;
centers = reshape(settings.channelCentersNm,1,[]);
sigma = reshape(settings.channelSigmaNm,1,[]);
cfg.channelResponse = exp(-0.5*((wavelengthNm(:)-centers)./sigma).^2);
periodsPerAxis = [4, 8, 8];
diameterMm = [0.20, 0.20, 0.40];
for k = 1:numel(periodsPerAxis)
    cfg.periods = repmat(periodsPerAxis(k),1,2);
    cfg.pupilRadius = diameterMm(k)*0.5e-3;
    fprintf('Trial %d: %d x %d pitches, pupil %.2f mm\n', k,cfg.periods,diameterMm(k));
    [optics,images] = runUDCImagingV2(cfg,stack,imageDir,scene,"srgb");
    x = optics.coordinate;
    frequency = optics.spatialFrequency;
    psf = optics.psf;
    mtf = optics.mtf;
    fixedExposure = images.fixedExposure;
    captured = optics.channelCapturedThroughput;
    openSamples = nnz(optics.pupil & optics.layerGeometry.unionMask);
    pupilSamples = nnz(optics.pupil);
    sampleArea = optics.layerGeometry.gridStep^2;
    incidentArea = pupilSamples*sampleArea;
    save(fullfile(out,sprintf('trial%d.mat',k)), 'cfg','stack','x','frequency','psf','mtf', ...
        'fixedExposure','captured','openSamples','pupilSamples','incidentArea','-v7');
end
