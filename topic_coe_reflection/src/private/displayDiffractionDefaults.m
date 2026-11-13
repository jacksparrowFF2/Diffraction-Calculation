function cfg = displayDiffractionDefaults(mode)
%displayDiffractionDefaults Parameters for COE reflection or UDC imaging (SI).
% Geometry follows CustomAperture_fraunhofer_fft: 50 um pitch, 10 um radius.
% All R/T values below are illustrative INTENSITY coefficients, not measurements.
    arguments
        mode (1,1) string {mustBeMember(mode, ["reflection", "udc"])} = "reflection"
    end
    cfg.mode = mode;
    cfg.pitch = 50e-6;
    cfg.radius = 10e-6;
    cfg.aspectRatio = 1;
    cfg.cellAnglesDeg = [30, 0; 0, -30];
    cfg.periods = [16, 16]; % [Nx, Ny], both multiples of two
    cfg.samplesPerPeriod = 32;
    cfg.paddingFactor = 4;
    cfg.propagationDistance = 10; % Reflection receiver distance z / m (Fraunhofer).
    cfg.incidentIrradiance = 1; % Per-wavelength reference incident irradiance / W m^-2.
    cfg.receiverBinSamples = 3; % Minimum odd midpoint samples per receiver pixel axis.
    cfg.source = coeSourceDefaults();
    cfg.mask = []; % Optional Ny-by-Nx sampled aperture, replaces customPattern.
    cfg.heightMap = 0; % Scalar or same size as mask; physical height in metres.
    cfg.deltaIndex = 0.20;
    cfg.wavelengths = [620, 550, 450]*1e-9;
    cfg.spectralWeights = [1, 1, 1]; % Include source SPD AND integration-bin widths.
    cfg.channelResponse = eye(3); % wavelength-by-channel; narrowband RGB proxies
    cfg.featureIntensity = 0.60;
    cfg.backgroundIntensity = 0.01;
    cfg.complexCoefficient = []; % Optional r or t: scalar, Ny-by-Nx, or Ny-by-Nx-by-L.
    cfg.pupilRadius = Inf;
    cfg.focalLength = 3e-3; % Used ONLY by the camera branch, not as propagation distance.
    cfg.sensorPitch = 1e-6;
    cfg.outputSize = 401; % Odd: unambiguous centre for convolution and PSF.
    cfg.angleHalfRangeDeg = 6; % Reflection uses paraxial angle coordinate.
    cfg.metricAngleDeg = 2.5; % Cone HALF-angle; 2.5 deg corresponds to a 5 deg full cone.
    cfg.zeroOrderRadiusFactor = 0.45;
    if mode == "udc"
        cfg.featureIntensity = 0.65;
        cfg.backgroundIntensity = 0.002;
        cfg.pupilRadius = 0.35e-3;
    end
end
