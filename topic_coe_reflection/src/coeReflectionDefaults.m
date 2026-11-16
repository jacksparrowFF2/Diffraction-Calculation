function cfg = coeReflectionDefaults()
%coeReflectionDefaults Create a complete COE reflection configuration in SI.
    arguments
    end
    % 反射默认值集中定义；R 为强度比例，材料常数为占位示例。
    cfg.mode = "reflection";
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
    cfg.source = coreSourceDefaults();
    cfg.mask = []; % Optional Ny-by-Nx sampled aperture, replaces coreCustomPattern.
    cfg.heightMap = 0; % Scalar or same size as mask; physical height in metres.
    cfg.deltaIndex = 0.20;
    cfg.wavelengths = [620, 550, 450]*1e-9;
    cfg.spectralWeights = [1, 1, 1]; % Include source SPD AND integration-bin widths.
    cfg.channelResponse = eye(3); % wavelength-by-channel; narrowband RGB proxies
    cfg.featureIntensity = 0.60;
    cfg.backgroundIntensity = 0.01;
    cfg.complexCoefficient = []; % Optional r: scalar, Ny-by-Nx, or Ny-by-Nx-by-L.
    cfg.pupilRadius = Inf;
    cfg.focalLength = 3e-3; % Retained for historical configuration compatibility.
    cfg.sensorPitch = 1e-6;
    cfg.outputSize = 1201; % Odd: unambiguous centre for convolution and PSF.
    cfg.angleHalfRangeDeg = 6; % Reflection uses paraxial angle coordinate.
    cfg.metricAngleDeg = 2.5; % Cone HALF-angle; 2.5 deg corresponds to a 5 deg full cone.
    cfg.zeroOrderRadiusFactor = 0.45;
end
