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
    cfg.periods = [4, 4]; % [Nx, Ny], both multiples of two
    cfg.samplesPerPeriod = 32;
    cfg.paddingFactor = 4;
    cfg.propagationDistance = 0.25; % Reflection receiver distance z / m (Fraunhofer).
    cfg.incidentIrradiance = 1; % Per-wavelength reference incident irradiance / W m^-2.
    cfg.receiverBinSamples = 3; % Minimum odd midpoint samples per receiver pixel axis.
    cfg.receiverView.enabled = false; % 复用远场角图，输出多个光栅→接收面距离的毫米坐标。
    cfg.receiverView.distances = [0.25]; % m；仅远场坐标/辐照度换算，不是近场重算。
    cfg.source = coreSourceDefaults();
    % 可选几何观察视图；与propagationDistance对应的远场接收面分开。
    cfg.observer.enabled = false;
    cfg.observer.model = "geometric"; % geometric / coherentEye：逐源点完整标量成像链。
    cfg.observer.distance = 0.25; % 眼睛到屏幕的法向距离 b / m。
    cfg.observer.pupilDiameter = 3e-3; % 人眼入口瞳直径；与屏幕pupilRadius独立。
    cfg.observer.pupilModel = "geometryOnly"; % geometryOnly / airyIntensity（非相干角图近似）
    cfg.observer.airyRadiusFactor = 20; % PSF截断半径/第一零点半径；不补回核外功率。
    cfg.observer.referenceMethod = "quadrature"; % quadrature / fftBins（宽面源比较）
    cfg.observer.referencePaddingFactor = 1; % FFT能量分箱的零填充倍数，正整数。
    cfg.observer.angleHalfRangeDeg = 2.5; % 固定的视野半角。
    cfg.observer.runBenchmark = false;
    cfg.observer.scanDistances = [0.15,0.20,0.25,0.30,0.35]; % m。
    cfg.observer.focusMode = "virtualSource"; % virtualSource=a+b / screen=b / distance。
    cfg.observer.focusDistance = 0.5; % focusMode=distance时生效，m。
    cfg.observer.retinaDistance = 17e-3; % 理想薄透镜到像面的距离，m；不代表眼睛处方。
    cfg.observer.pupilSamples = 1537; % 瞳面每轴点数，奇数；需满足复场/输出角域采样。
    cfg.observer.imageZoomHalfAngleDeg = 0.25;
    cfg.observer.displayReferenceHalfAngleDeg = 0.25; % 仅成像显示的共同曝光参考ROI；0取完整图。
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
