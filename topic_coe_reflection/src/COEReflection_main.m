%% COE 反射衍射：输入 → 参数 → 计算 → 展示 → 保存
% 直接运行本脚本；常用参数集中在前两节，具体 FFT 和能量积分算法在 private/ 中。
% 模型：法向局部相干膜堆 + 标量近轴远场；额外 heightMap 默认 0 nm。
% 默认 Diamond 图样和 BM/CF 的 n/k、厚度为可调示例，不代表实际产品。光谱采用 CIE D65。

%% 0. 初始化运行环境
clearvars;
clc;
close all;
% 仅清理工作区、命令窗口和旧图窗，保留现有搜索路径及 MCP 连接。

%% 1. 文件输入与输出设置
moduleDir = fileparts(fileparts(mfilename('fullpath')));
inputRoot = fullfile(moduleDir, 'Input');
projectName = "Case1";          % 项目文件夹名称：Case1=PNG，Case2=BMP，Case3=JPG；空字符串使用 Input 根目录
assert(~any(projectName == [".", ".."]) && ~contains(projectName, ["/", "\", ":"]), ...
    'COEReflection:ProjectName', '项目名称必须是单个文件夹名称，不能包含路径分隔符。');
useProjectConfig = true;       % 读取 Input/项目名称/config.txt；false 使用下方脚本参数
inputDir = fullfile(inputRoot, projectName);
outputDir = fullfile(moduleDir, 'Output', projectName); % 各项目结果分开保存
inputFile = "";                  % 旧单层/额外高度/外部复系数 MAT；实际多层图样请用 stack.geometry
caseName = "coe_reflection";     % 输出文件名前缀；同名结果会覆盖
saveResults = true;             % 是否保存 MAT、CSV 和 PNG

% MAT 文件可包含以下任意字段，读取后覆盖脚本中的对应值：
% mask：二维 0~1 材料区域图，推荐二值；1 为特征区，0 为背景区。
% heightMap：高度（m），标量或与 mask 同尺寸的矩阵。
% complexCoefficient：复振幅反射系数 r，标量、二维矩阵或 Ny×Nx×波长数数组。
% 这里输入的是反射面数据，照明为平面波，不需要场景图片。

%% 2. 模型参数（长度统一使用 m，角度标有 Deg 的使用 °）
cfg = coeReflectionDefaults();   % 补齐底层所需字段；本脚本列出反射计算的常用参数
opticalModel = "multilayer";     % multilayer：阳极/BM/CF；simple：保留旧单界面模型
mustBeMember(opticalModel, ["multilayer", "simple"]);

% 2.1 几何与采样
cfg.pitch = 50e-6;              % 像素周期，50 μm
cfg.periods = [16, 16];         % [x 方向周期数, y 方向周期数]，均为偶数
cfg.samplesPerPeriod = 64;      % 每周期采样点数；采样间距 = pitch / samplesPerPeriod
cfg.paddingFactor = 2;          % 兼容旧配置；当前反射传播的输出采样由 outputSize/角域控制
cfg.propagationDistance = 10;   % 接收面距离 z，单位 m；须满足夫琅禾费远场条件
cfg.incidentIrradiance = 1;     % 每个波长单独以此入射辐照度标定，单位 W/m^2
cfg.receiverBinSamples = 3;     % 每像素每轴最少积分点数（奇数），窄峰时自动增加
cfg.pupilRadius = Inf;          % Inf 表示整个矩形计算面参与反射
cfg.mask = [];                 % 留空时自动生成几何；自定义矩阵尺寸见下方 gridSize

% 2.2 光谱与反射面
spectrumMode = "D65";          % "D65"：380~780 nm 色度合成；"RGB"：原三波长演示
wavelengthStepNm = 5;          % 光谱步长；可改为 1/2/5/10 nm，越细计算越久
previewWavelengthNm = [620, 550, 450]; % 远场图仅预览这些波长，全部波长均参与积分
cfg.wavelengths = [620, 550, 450]*1e-9; % RGB 模式下的三个代理波长
cfg.spectralWeights = [1, 1, 1]; % RGB 模式使用的入射权重
cfg.channelResponse = eye(3);   % RGB 模式的 R/G/B 通道
% 以下为旧单界面模型参数，仅 opticalModel="simple" 时使用。
cfg.materials.enabled = true;  % true：由 n/k 算复反射系数；false：用下面旧的强度反射率
cfg.materials.incidentN = 1.0;  % 入射介质为空气
cfg.materials.featureN = 1.8;   % 特征区材料示例 n，整个波段暂为常数
cfg.materials.featureK = 0.05;  % 特征区消光系数 k
cfg.materials.backgroundN = 1.5;
cfg.materials.backgroundK = 0.0;
% 法向单界面 r=(n_incident-(n+i*k))/(n_incident+(n+i*k))，R=abs(r)^2。
% 外部 complexCoefficient 非空时优先于材料 n/k，适合以后输入薄膜计算结果。
cfg.featureIntensity = 0.60;    % 特征区强度反射率 R；算法使用 sqrt(R) 作为振幅
cfg.backgroundIntensity = 0.01; % 背景区强度反射率；也可按波长指定一组值
cfg.deltaIndex = 0.20;          % 高度相位对应的折射率差 Δn
cfg.heightMap = 0;             % 默认额外高度 0 nm；暂不添加高度相位
cfg.complexCoefficient = [];   % 留空时由 mask 和反射率构造；非空时优先使用给定的 r
% 反射高度相位：phi = 4*pi*deltaIndex*heightMap/lambda。
% 均匀高度只改变整体相位，不改变远场强度；非均匀高度才会重新分配衍射能量。
% 若输入 r 已包含该高度相位，请将 heightMap 设为 0，避免重复计入。
% deltaIndex 仍是额外光程层的有效折射率差，不自动等同于任一种材料的 n。

% 2.2.1 图案化多层 COE：图样输入 + BM/CF 材料与厚度
stack = coeStackDefaults();
stack.geometry.source = "images"; % images：项目图片集（默认）；mat：同一 MAT
stack.geometry.matFile = "coe_layers.mat"; % 字段 anodeMask、cfR/G/BMask，可选 bmOpenMask 或 bmMask
stack.geometry.anodeFile = "anode"; % 不写扩展名时自动匹配 PNG/BMP/JPG/JPEG/TIF/TIFF
stack.geometry.bmFile = "BM";  % 空字符串时自动从阳极外扩生成 BM 开口
stack.geometry.bmImageIsOpening = false; % 示例 BM 图片白色=实心材料；白色=开口时改 true
stack.geometry.cfFiles = ["CF_R", "CF_G", "CF_B"];
stack.geometry.threshold = 0.5; % 图片二值化阈值；推荐同尺寸、同原点的灰度黑白 PNG，无透明通道
stack.geometry.requireCFCoverage = true; % CF 必须覆盖 BM 开口，各色之间不可交叠
% 实测输入：改 source="mat" 或 "images"，文件名均相对于 Input/projectName/，可含子目录。
% 输入必须覆盖完整计算面：[periods(2),periods(1)]*samplesPerPeriod，行 y、列 x。
% 不自动缩放/平铺；输入 BM 后直接采用该图样，bmExpansion 只在缺省 BM 时生成开口。
% 详细字段、物理尺寸、nk CSV 示例见 docs/README.md。
stack.bmExpansion = 5e-6;       % BM 开口相对阳极边界向外扩张 5 μm，可调
stack.anodeReflectance = 1.0;   % 仅阳极图样内的裸底为 100% 反射
stack.anodePhaseDeg = 180;      % 理想电反射边界相位；阳极外为匹配吸收底（裸底 r=0）
stack.materials.BM.n = 1.8; stack.materials.BM.k = 0.35; stack.materials.BM.thickness = 1.0e-6;
stack.materials.R.n = 1.65; stack.materials.R.k = 0.03; stack.materials.R.thickness = 1.1e-6;
stack.materials.G.n = 1.60; stack.materials.G.k = 0.04; stack.materials.G.thickness = 1.0e-6;
stack.materials.B.n = 1.55; stack.materials.B.k = 0.05; stack.materials.B.thickness = 0.9e-6;
stack.materials.BM.nkFile = ""; % 实测时填写如 "Materials/BM_nk.csv"，空字符串使用上面的 n/k
stack.materials.R.nkFile = "";  % 如 "Materials/CF_R_nk.csv"
stack.materials.G.nkFile = "";  % 如 "Materials/CF_G_nk.csv"
stack.materials.B.nkFile = "";  % 如 "Materials/CF_B_nk.csv"
% 每种材料还可设 nkFile="BM_nk.csv" 等：三列 [波长 nm, n, k]，文件优先于上面的常数。
% 也可将 n/k 设为与 cfg.wavelengths 一一对应的向量。常数占位值不保证 CF 真正具有 RGB 通带。
% CF 在 BM 上方的部分为 CF/BM/吸收底；BM 外扩开口环带为 CF/吸收底，绝不扩大阳极反射范围。
% 膜堆厚度和共同参考平面传播相位属于多层模型本身，两种比较状态均保留。

% 2.2.2 非均匀高度示例（仅 inputFile 为空时启用；实际文件输入优先）
useHeightExample = false;       % 默认关闭；需要演示非均匀高度时再改为 true
heightExample.height = 687.5e-9; % 条纹高区 687.5 nm，低区为 0；550 nm 光下相位差为 pi
heightExample.period = cfg.pitch/4; % 高度条纹周期 12.5 μm，与材料区域的 50 μm 周期不同
heightExample.dutyCycle = 0.5;  % 每周期高区占比 50%
heightExample.direction = "x"; % "x"：高度沿 x 交替，条纹平行 y；"y" 则相反
% 示例作用于整个计算面（特征区和背景区），材料几何及强度反射率保持不变。
% 这是验证相位作用的演示结构，不代表实际产品或优化结果。
% 预计一级偏转约 lambda/period：RGB 分别为 2.84°、2.52°、2.06°。

% 2.3 观察窗口、评价范围与显示
cfg.outputSize = 401;           % 输出角度网格的边长，必须为奇数
cfg.angleHalfRangeDeg = 6;      % 观察角度范围 [-6°, 6°]，近轴模型上限为 10°
cfg.metricAngleDeg = 2.5;       % 评价圆锥半角 2.5°，对应全角 5°，中心为镜面反射方向
cfg.zeroOrderRadiusFactor = 0.45; % 零级积分半径 = 此系数 * lambda / pitch，单位 rad
% 独立于雾度角度的观察角域；中心 [theta_x, theta_y]，圆锥/环带必须完整位于输出窗口内。
observation.centerDeg = [0, 0];
observation.halfAngleDeg = 2.5; % 外半角，可按实际探测器或人眼接收视角调整
observation.innerHalfAngleDeg = 0; % 0=完整圆锥；>0=排除中心的环带
observation.curveSamples = 301; % 累计曲线的半径采样数
cfg.observation = observation; % 保存到本次配置，供图中标注
% 非零级衍射斑：先从原始能量图共同定位，再对两个状态使用同一积分圆比较。
spotOptions.zeroOrderExclusionDeg = 0.30; % 以镜面方向 (0,0) 为中心的零级排除半角
spotOptions.integrationHalfAngleDeg = 0.15; % 每个光斑的圆形积分半角（非直径）
spotOptions.minSeparationDeg = 0.35; % 光斑最小中心间距；至少为两倍积分半角以避免重叠
spotOptions.maxSpots = 8;       % 每个波长 / 宽带组最多保留的主要光斑数
spotOptions.minRelativePeak = 0.08; % 相对零级排除区外最强峰的检测阈值
spotOptions.wavelengthNm = previewWavelengthNm; % [] 可统计全部波长；另始终统计宽带辐射度/明视觉
displayFloorDb = -60;           % 图像显示下限；仅影响绘图，不改变物理数组和指标
displayStyle.fontName = 'Microsoft YaHei'; % 中文、英文和数字使用同一字体
displayStyle.fontSize = 11;
displayStyle.lineWidth = 1;
% 图形标记开关：只控制绘图，不影响检测、积分、CSV 或原始 XYZ。
displayStyle.showAnnotations = true; % 总开关：false 隐藏下列所有辅助标记
displayStyle.showMetricCone = true; % 青色：雾度评价圆锥边界
displayStyle.showObservationRegion = true; % 黄色：观察角域；也控制累计曲线上的半角参考线
displayStyle.showZeroOrderExclusion = true; % 红色：零级排除圆
displayStyle.showSpotRegions = true; % 品红色：非零级光斑积分圆
displayStyle.showSpotLabels = true; % S 编号和对应引线
rgbOptions.clipQuantile = 0.9985; % 参考孔径模块：亮部截断 + 对数压缩 + 伽马增强
rgbOptions.alpha = 60;
rgbOptions.gamma = 0.72;
colorOptions.fixedGain = 0;    % 0=两状态共同自动曝光；>0=固定 XYZ 密度→线性 RGB 的曝光倍数
colorOptions.exposureEV = 0;   % 在共同曝光基础上的 EV 调整；+1 表示线性亮度乘 2
colorOptions.toneMapping = "log"; % "log"：增强弱衍射；"none"：关闭增强，恢复线性曝光图
colorOptions.toneStrength = 3000; % 对数亮度增强强度，越大越容易看见弱衍射
colorOptions.luminanceGamma = 0.75; % 仅作用于 Y；小于 1 提亮暗部，1 为不额外加亮度伽马
colorView.centerDeg = [0, 0];  % D65 专用图下排局部放大的中心 [theta_x, theta_y]
colorView.halfAngleDeg = 1.5;  % 下排显示 ±1.5°，仅裁切显示，不修改积分角域或重新计算
phaseStateLabels = ["不叠加高度相位", "叠加高度相位"];

%% 2.4 项目配置（覆盖上面的默认参数；在光谱/几何计算之前生效）
projectConfig = struct();
if useProjectConfig
    settings = struct();
    settings.cfg = cfg;
    settings.stack = stack;
    settings.opticalModel = opticalModel;
    settings.spectrumMode = spectrumMode;
    settings.wavelengthStepNm = wavelengthStepNm;
    settings.previewWavelengthNm = previewWavelengthNm;
    settings.inputFile = inputFile;
    settings.caseName = caseName;
    settings.saveResults = saveResults;
    settings.useHeightExample = useHeightExample;
    settings.heightExample = heightExample;
    settings.observation = observation;
    settings.spotOptions = spotOptions;
    settings.displayFloorDb = displayFloorDb;
    settings.displayStyle = displayStyle;
    settings.rgbOptions = rgbOptions;
    settings.colorOptions = colorOptions;
    settings.colorView = colorView;
    [settings, projectConfig] = loadCoeProjectConfig(fullfile(inputDir, 'config.txt'), settings);
    cfg = settings.cfg;
    stack = settings.stack;
    opticalModel = settings.opticalModel;
    spectrumMode = settings.spectrumMode;
    wavelengthStepNm = settings.wavelengthStepNm;
    previewWavelengthNm = settings.previewWavelengthNm;
    inputFile = settings.inputFile;
    caseName = settings.caseName;
    saveResults = settings.saveResults;
    useHeightExample = settings.useHeightExample;
    heightExample = settings.heightExample;
    observation = settings.observation;
    spotOptions = settings.spotOptions;
    displayFloorDb = settings.displayFloorDb;
    displayStyle = settings.displayStyle;
    rgbOptions = settings.rgbOptions;
    colorOptions = settings.colorOptions;
    colorView = settings.colorView;
    fprintf('读取项目配置：%s（%d 项覆盖）。\n', projectConfig.file, numel(projectConfig.keys));
end
cfg.observation = observation;
mustBeMember(opticalModel, ["multilayer", "simple"]);
mustBeMember(spectrumMode, ["D65", "RGB"]);
if spectrumMode == "D65"
    validateattributes(wavelengthStepNm, {'numeric'}, {'scalar', 'integer', 'positive', '<=', 20});
    assert(mod(400, wavelengthStepNm) == 0, 'COEReflection:SpectrumStep', '波长步长必须整除 400 nm。');
    spectrum = loadCoeSpectrum(380:wavelengthStepNm:780);
    cfg.wavelengths = spectrum.wavelengthNm.'/1e9;
    cfg.spectralWeights = spectrum.energyWeights.'; % D65 SPD × 梯形积分步长，归一化入射能量
    cfg.channelResponse = ones(numel(cfg.wavelengths), 1); % 这里只形成宽带能量；色度在第 6 节单独积分
    spectrumLabel = "D65 380–780 nm";
    photopicLabel = "D65 明视觉";
else
    spectrum.energyWeights = cfg.spectralWeights(:)/sum(cfg.spectralWeights);
    spectrum.photopicWeights = spectrum.energyWeights; % 旧模式不作标准光度计算
    spectrumLabel = "三波长等权代理";
    photopicLabel = "三波长代理（非标准光度）";
end

%% 3. 读取设计文件并检查输入尺寸
gridSize = [cfg.periods(2), cfg.periods(1)]*cfg.samplesPerPeriod;
if strlength(inputFile) > 0
    inputPath = fullfile(inputDir, inputFile);
    assert(isfile(inputPath), 'COEReflection:MissingInput', '输入文件不存在：%s', inputPath);
    inputData = load(inputPath, 'mask', 'heightMap', 'complexCoefficient');
    assert(~isempty(fieldnames(inputData)), 'COEReflection:EmptyInput', ...
        'MAT 文件应至少包含 mask、heightMap 或 complexCoefficient 之一。');
    if isfield(inputData, 'mask')
        cfg.mask = inputData.mask;
    end
    if isfield(inputData, 'heightMap')
        cfg.heightMap = inputData.heightMap;
    end
    if isfield(inputData, 'complexCoefficient')
        cfg.complexCoefficient = inputData.complexCoefficient;
    end
    fprintf('读取设计文件：%s\n', inputPath);
else
    fprintf('未使用额外 MAT 覆盖，按主脚本的图样与材料设置计算。\n');
    if useHeightExample
        cfg.heightMap = createCoeHeightExample(cfg, heightExample); %#ok<UNRCH> 可选示例默认关闭
        fprintf('二值高度条纹：高区 %.1f nm，周期 %.2f μm，占空比 %.2f，变化方向 %s。\n', ...
            heightExample.height*1e9, heightExample.period*1e6, ...
            heightExample.dutyCycle, heightExample.direction);
        fprintf('各波长高低区相位差 / pi：');
        fprintf(' %.3f', 4*cfg.deltaIndex*heightExample.height./cfg.wavelengths);
        fprintf('\n');
    end
end
% 所有二维输入覆盖完整计算面，行为 y，列为 x；不自动缩放或插值。
assert(isempty(cfg.mask) || isequal(size(cfg.mask), gridSize), ...
    'COEReflection:MaskSize', 'mask 尺寸应为 %d×%d（行×列）。', gridSize(1), gridSize(2));
assert(isscalar(cfg.heightMap) || isequal(size(cfg.heightMap), gridSize), ...
    'COEReflection:HeightSize', 'heightMap 应为标量或与 mask 同尺寸的矩阵。');
fprintf('计算面采样：%d×%d，采样间距 %.3f μm。\n', gridSize(1), gridSize(2), ...
    cfg.pitch/cfg.samplesPerPeriod*1e6);

%% 3.1 建立局部膜层区域与复反射系数（不依赖其它 topic）
stackGeometry = struct();
if opticalModel == "multilayer"
    assert(isempty(cfg.mask), 'COEReflection:LegacyMask', ...
        '多层模式请通过 stack.geometry 导入各层图样，不使用旧单一 mask。');
    [cfg.coeStack, stackGeometry] = prepareCoeStack(cfg, stack, inputDir);
    fprintf('多层图样来源：%s；项目目录：%s。\n', stack.geometry.source, inputDir);
    cfg.mask = double(stackGeometry.anodeMask); % FFT 坐标沿用现有接口；各层系数由区域表独立决定
    fprintf('多层 COE：%d 种局部膜堆，共同参考面高度 %.3f μm。\n', ...
        height(cfg.coeStack.regionTable), cfg.coeStack.referenceHeight*1e6);
    fprintf('阳极面积 %.2f%%，BM 外扩吸收环带 %.2f%%；额外 heightMap 保持用户设定值。\n', ...
        100*mean(stackGeometry.anodeMask, 'all'), 100*mean(stackGeometry.expansionAbsorberMask, 'all'));
    disp(cfg.coeStack.regionTable);
end

%% 4. 计算“叠加高度相位”状态的反射衍射
% 底层函数完成：几何构建 → 复振幅反射面 → 双程相位 → FFT → 角度网格能量积分。
% 不同波长按强度非相干叠加，kernel 保留反射损耗和有限观察窗口的能量损失。
result = simulateDisplayDiffraction(cfg);

%% 5. 计算“不叠加高度相位”状态并比较指标
% 这是同一输入结构的两种计算状态，第二种状态仅将额外 heightMap 置零。
% 两种状态均保留输入复系数自身的相位，因此不能称为“无相位 / 有相位”。
% baseline 变量及字段名沿用函数接口命名，图表使用上方明确的状态名称。
baselineCfg = cfg;
baselineCfg.heightMap = 0;
baseline = simulateDisplayDiffraction(baselineCfg);
result.baselineKernel = baseline.kernel;
result.baselineMetrics = baseline.metrics;
% 零级保持率 = 叠加高度相位 / 不叠加高度相位的零级能量占比；分母为零时记为 NaN。
retention = NaN(height(result.metrics), 1);
active = baseline.metrics.ZeroOrderFraction > 0;
retention(active) = result.metrics.ZeroOrderFraction(active)./baseline.metrics.ZeroOrderFraction(active);
result.metrics.ZeroOrderRetention = retention;
disp('叠加高度相位：逐波长指标');
disp(result.metrics(1:max(1, floor(numel(cfg.wavelengths)/8)):end, :));
disp('不叠加高度相位：逐波长指标');
disp(baseline.metrics(1:max(1, floor(numel(cfg.wavelengths)/8)):end, :));
% Throughput：总反射能量 / 入射能量；CapturedFraction：窗口内能量 / 总反射能量。
% ZeroOrderFraction：零级能量 / 总反射能量。
% OutsideAngleFraction：评价圆锥外能量 / 总反射能量，包含窗口外能量。
% 用其百分比作为反射雾度评价值；这是模型代理指标，不是 ASTM D1003 透射雾度。
% ParsevalAbsoluteError：FFT 能量守恒的绝对误差。
% 关闭高度示例且 heightMap=0 时两种状态相同，零级保持率应为 1（有非零分母时）。
% 明确记录状态、半角、全角与百分比，避免将 5° 全角误用为 5° 半角。
rowCount = numel(cfg.wavelengths);
comparisonMetrics = table( ...
    [repmat(phaseStateLabels(1), rowCount, 1); repmat(phaseStateLabels(2), rowCount, 1)], ...
    repmat(cfg.wavelengths(:)*1e9, 2, 1), repmat(cfg.metricAngleDeg, 2*rowCount, 1), ...
    repmat(2*cfg.metricAngleDeg, 2*rowCount, 1), ...
    100*[baseline.metrics.OutsideAngleFraction; result.metrics.OutsideAngleFraction], ...
    100*[baseline.metrics.ZeroOrderFraction; result.metrics.ZeroOrderFraction], ...
    [baseline.metrics.CapturedFraction; result.metrics.CapturedFraction], ...
    'VariableNames', {'HeightPhaseState', 'WavelengthNm', 'ConeHalfAngleDeg', ...
    'ConeFullAngleDeg', 'ReflectionHazeProxyPercent', 'ZeroOrderPercent', 'CapturedFraction'});
% 同时给出以入射能量为分母的绝对反射率，区分“反射变少”与“反射能量被移到其它方向”。
comparisonMetrics.TotalReflectancePercent = 100*[baseline.metrics.Throughput; result.metrics.Throughput];
comparisonMetrics.ConeOutsideReflectancePercent = comparisonMetrics.TotalReflectancePercent ...
    .*comparisonMetrics.ReflectionHazeProxyPercent/100;
comparisonMetrics.ZeroOrderReflectancePercent = comparisonMetrics.TotalReflectancePercent ...
    .*comparisonMetrics.ZeroOrderPercent/100;
disp('高度相位状态对比（雾度评价值为模型反射锥外占比，不代表标准仪器雾度）：');
disp(comparisonMetrics(1:max(1, floor(rowCount/8)):end, :));

%% 5.1 观察角域绝对反射率与累计角能量
% 每个波长的 spectralEnergy 已除以入射能量；先积分能量，再作光谱加权，不能平均各波长占比。
withoutAngular = computeCoeAngularMetrics(baseline, observation, spectrum);
withAngular = computeCoeAngularMetrics(result, observation, spectrum);
comparisonMetrics.ViewAbsoluteReflectancePercent = 100*[withoutAngular.spectralViewReflectance; ...
    withAngular.spectralViewReflectance];
observationMetrics = table(phaseStateLabels.', ...
    100*[withoutAngular.radiantViewReflectance; withAngular.radiantViewReflectance], ...
    100*[withoutAngular.photopicViewReflectance; withAngular.photopicViewReflectance], ...
    100*[withoutAngular.radiantTotalReflectance; withAngular.radiantTotalReflectance], ...
    100*[withoutAngular.photopicTotalReflectance; withAngular.photopicTotalReflectance], ...
    'VariableNames', {'HeightPhaseState', 'ViewRadiantReflectancePercent', ...
    'ViewPhotopicReflectancePercent', 'TotalRadiantReflectancePercent', 'TotalPhotopicReflectancePercent'});
observationMetrics.CenterXDeg = repmat(observation.centerDeg(1), 2, 1);
observationMetrics.CenterYDeg = repmat(observation.centerDeg(2), 2, 1);
observationMetrics.InnerHalfAngleDeg = repmat(observation.innerHalfAngleDeg, 2, 1);
observationMetrics.OuterHalfAngleDeg = repmat(observation.halfAngleDeg, 2, 1);
if spectrumMode == "D65"
    viewXYZ = [withoutAngular.spectralViewReflectance.'; withAngular.spectralViewReflectance.']*spectrum.xyzWeights;
    totalXYZ = sum(viewXYZ, 2);
    totalXYZ(totalXYZ == 0) = NaN;
    observationMetrics.ViewChromaticityX = viewXYZ(:, 1)./totalXYZ;
    observationMetrics.ViewChromaticityY = viewXYZ(:, 2)./totalXYZ;
end
cumulativeMetrics = table(withAngular.radiusDeg, ...
    100*withoutAngular.radiantCumulativeReflectance, 100*withAngular.radiantCumulativeReflectance, ...
    100*withoutAngular.photopicCumulativeReflectance, 100*withAngular.photopicCumulativeReflectance, ...
    'VariableNames', {'RadiusDeg', 'WithoutHeightRadiantPercent', 'WithHeightRadiantPercent', ...
    'WithoutHeightPhotopicPercent', 'WithHeightPhotopicPercent'});
disp('观察角域内绝对反射率（分母为入射能量 / 入射明视觉光通量，单位 %）：');
disp(observationMetrics);

%% 5.2 非零级光斑位置与绝对反射率
% S1、S2 等仅是检测到的局部亮斑编号，不自动代表理论 (m,n) 衍射阶数。
% 全部统计基于未增强的物理数组，分母为入射能量（photopic 组为入射明视觉光通量）。
diffractionSpots = computeCoeDiffractionSpots(baseline, result, spectrum, spotOptions);
spotMetrics = diffractionSpots.metrics;
disp('非零级主要衍射斑（宽带明视觉组；RGB 模式为代理权重）：');
disp(spotMetrics(spotMetrics.Group == "photopic", :));

%% 6. 展示反射面与各波长的远场结果
geometryFigure = figure('Name', 'COE 反射面', 'Color', 'w', 'Position', [80, 80, 1050, 430]);
if opticalModel == "multilayer"
    set(geometryFigure, 'Position', [60, 80, 1450, 430]);
    tiledlayout(1, 4, 'TileSpacing', 'compact');
else
    tiledlayout(1, 2);
end
nexttile;
imagesc(result.x*1e6, result.y*1e6, result.mask);
axis image xy;
colorbar;
clim([0, 1]);
xlabel('x / μm'); ylabel('y / μm'); title('阳极 / 旧模型特征区');
if opticalModel == "multilayer"
    nexttile;
    imagesc(result.x*1e6, result.y*1e6, stackGeometry.bmMask);
    axis image xy; colorbar; clim([0, 1]);
    xlabel('x / μm'); ylabel('y / μm'); title('BM 实体（开口内含吸收环带）');
    nexttile;
    imagesc(result.x*1e6, result.y*1e6, stackGeometry.cfLabel);
    axis image xy; clim([-0.5, 3.5]);
    colormap(gca, [0.1, 0.1, 0.1; 1, 0.2, 0.2; 0.2, 0.8, 0.2; 0.2, 0.4, 1]);
    cb = colorbar; cb.Ticks = 0:3; cb.TickLabels = {'无 CF', 'R', 'G', 'B'};
    xlabel('x / μm'); ylabel('y / μm'); title('CF 图样（标签色，非光学色）');
end
nexttile;
heightNm = cfg.heightMap.*ones(size(result.mask))*1e9;
imagesc(result.x*1e6, result.y*1e6, heightNm);
axis image xy;
colorbar;
xlabel('x / μm'); ylabel('y / μm'); title('额外高度 / nm');

% 横向排列：每列一个波长 / RGB；上行为不叠加、下行为叠加高度相位。
% 单波长图共用能量参考和灰度 dB 色标；RGB 使用独立的共同显示映射，无 dB 色标。
angleDeg = rad2deg(result.coordinate);
referenceEnergy = max([baseline.spectralEnergy(:); result.spectralEnergy(:)]);
referenceEnergy = max(referenceEnergy, realmin);
hasRGB = spectrumMode == "D65" || size(cfg.channelResponse, 2) == 3;
displayData = struct('referenceEnergy', referenceEnergy, 'floorDb', displayFloorDb);
colorData = struct();
if spectrumMode == "D65"
    [colorData.withoutHeight, colorData.withHeight, colorData.mapping] = ...
        synthesizeCoeColor(baseline, result, spectrum, colorOptions);
    displayData.withoutHeightRGB = colorData.withoutHeight.displaySRGB;
    displayData.withHeightRGB = colorData.withHeight.displaySRGB;
    rgbTitle = "D65 / sRGB（显示增强）";
    if colorOptions.toneMapping == "none"
        rgbTitle = "D65 / sRGB（线性曝光）";
    end
    fprintf('CIE D65 白点 xy = %.5f, %.5f；共同曝光倍数 %.6g。\n', spectrum.whiteXY, colorData.mapping.gain);
    fprintf('两状态 sRGB 色域/高光裁剪像素比例：%.2f%% / %.2f%%。\n', ...
        100*colorData.withoutHeight.gamutClippedFraction, 100*colorData.withHeight.gamutClippedFraction);
elseif hasRGB
    % channelSignal 已计入光谱权重及响应；channelResponse 的列应依次表示 R/G/B。
    [displayData.withoutHeightRGB, displayData.withHeightRGB, displayData.rgbMapping] = ...
        composeCoeReflectionRgb(baseline.channelSignal, result.channelSignal, rgbOptions);
    rgbTitle = "三波长 RGB 代理（增强显示）";
else
    warning('COEReflection:RGBChannels', '通道数不是 3，跳过 RGB 预览；单波长图正常输出。');
end
diffractionFigure = figure('Name', 'COE 远场反射：高度相位作用', 'Color', 'w', ...
    'Units', 'normalized', 'Position', [0.03, 0.07, 0.94, 0.82]);
previewIndices = zeros(size(previewWavelengthNm));
for j = 1:numel(previewWavelengthNm)
    [~, previewIndices(j)] = min(abs(cfg.wavelengths*1e9 - previewWavelengthNm(j)));
end
previewIndices = unique(previewIndices, 'stable');
% 两种坐标图直接复用同一组 dB 数组：同参考、同下限、同波长、同高度状态。
% 等面积像素下 E/Eref=I/Iref；原始 W/m^2 数组仍保存在 result 中。
displayData.referenceIrradiance = referenceEnergy*result.incidentPower/result.receiverPixelArea;
displayData.previewWavelengthNm = cfg.wavelengths(previewIndices)*1e9;
displayData.withoutHeightDb = 10*log10(max( ...
    baseline.spectralEnergy(:, :, previewIndices)/referenceEnergy, 10^(displayFloorDb/10)));
displayData.withHeightDb = 10*log10(max( ...
    result.spectralEnergy(:, :, previewIndices)/referenceEnergy, 10^(displayFloorDb/10)));
stateDbImages = {displayData.withoutHeightDb, displayData.withHeightDb};
columnCount = numel(previewIndices) + double(hasRGB);
diffractionLayout = tiledlayout(2, columnCount, ...
    'TileSpacing', 'compact', 'Padding', 'compact');
title(diffractionLayout, '同一结构 · 高度相位作用对比');
annotationHint = "辅助标记已关闭";
if displayStyle.showAnnotations
    annotationHint = "辅助标记按 displayStyle 分类开关显示";
end
subtitle(diffractionLayout, spectrumLabel + "｜上/下：两种高度相位状态；" + annotationHint);
stateImages = {};
if hasRGB
    stateImages = {displayData.withoutHeightRGB, displayData.withHeightRGB};
end
for state = 1:2
    for column = 1:numel(previewIndices)
        k = previewIndices(column);
        nexttile((state-1)*columnCount + column);
        imagesc(angleDeg, angleDeg, stateDbImages{state}(:, :, column));
        formatFarFieldAxes(gca, cfg, displayStyle);
        colormap(gca, gray(256)); clim([displayFloorDb, 0]);
        cb = colorbar; cb.Ticks = linspace(displayFloorDb, 0, 4); cb.Label.String = '相对强度 / dB';
        title(sprintf('%s · %.0f nm', phaseStateLabels(state), cfg.wavelengths(k)*1e9));
        markDiffractionSpots(gca, diffractionSpots, "wavelength", cfg.wavelengths(k)*1e9, displayStyle);
    end
    if hasRGB
        nexttile(state*columnCount);
        imagesc(angleDeg, angleDeg, stateImages{state});
        formatFarFieldAxes(gca, cfg, displayStyle);
        title(phaseStateLabels(state) + " · " + rgbTitle);
        markDiffractionSpots(gca, diffractionSpots, "photopic", NaN, displayStyle);
    end
end

% 单独输出大尺寸 D65 图：上排保留全范围，下排放大指定区域。
colorFigure = gobjects(0);
if spectrumMode == "D65"
    validateattributes(colorView.centerDeg, {'numeric'}, {'size', [1, 2], 'real', 'finite'});
    validateattributes(colorView.halfAngleDeg, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
    assert(all(abs(colorView.centerDeg) + colorView.halfAngleDeg <= cfg.angleHalfRangeDeg), ...
        'COEReflection:ColorViewExtent', '色彩图放大范围必须位于输出角度窗口内。');
    colorFigure = figure('Name', 'D65 合成反射：全范围与局部放大', 'Color', 'w', ...
        'Position', [80, 60, 1120, 850]);
    colorLayout = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(colorLayout, 'D65 合成反射 · 两状态共用亮度映射');
    subtitle(colorLayout, sprintf('映射 %s｜强度 %g｜亮度 gamma %.2f｜曝光 EV %+.1f', ...
        colorOptions.toneMapping, colorOptions.toneStrength, colorOptions.luminanceGamma, colorOptions.exposureEV));
    colorImages = {displayData.withoutHeightRGB, displayData.withHeightRGB};
    for region = 1:2
        for state = 1:2
            nexttile;
            imagesc(angleDeg, angleDeg, colorImages{state});
            formatFarFieldAxes(gca, cfg, displayStyle);
            if region == 1
                title(phaseStateLabels(state) + " · 全范围");
            else
                xlim(colorView.centerDeg(1) + [-1, 1]*colorView.halfAngleDeg);
                ylim(colorView.centerDeg(2) + [-1, 1]*colorView.halfAngleDeg);
                xticks(colorView.centerDeg(1) + linspace(-1, 1, 5)*colorView.halfAngleDeg);
                yticks(colorView.centerDeg(2) + linspace(-1, 1, 5)*colorView.halfAngleDeg);
                title(phaseStateLabels(state) + " · 局部放大");
            end
            markDiffractionSpots(gca, diffractionSpots, "photopic", NaN, displayStyle);
        end
    end
end

metricsFigure = figure('Name', 'COE 高度相位作用与雾度评价', 'Color', 'w', ...
    'Position', [120, 120, 1150, 450]);
tiledlayout(1, 2);
% 使用序号作横坐标，保留输入顺序，也支持单波长或重复波长。
wavelengthAxis = cfg.wavelengths*1e9;
nexttile;
plot(wavelengthAxis, 100*baseline.metrics.ZeroOrderFraction, '-o', 'LineWidth', 1.5);
hold on;
plot(wavelengthAxis, 100*result.metrics.ZeroOrderFraction, '--x', 'LineWidth', 1.5);
xticks(sort(unique(cfg.wavelengths(previewIndices)*1e9))); xlabel('波长 / nm');
ylabel('零级能量占比 / %'); legend(phaseStateLabels, 'Location', 'best'); title('零级能量占比');
grid on; ylim([0, 100]);
nexttile;
plot(wavelengthAxis, 100*baseline.metrics.OutsideAngleFraction, '-o', 'LineWidth', 1.5);
hold on;
plot(wavelengthAxis, 100*result.metrics.OutsideAngleFraction, '--x', 'LineWidth', 1.5);
xticks(sort(unique(cfg.wavelengths(previewIndices)*1e9))); xlabel('波长 / nm');
ylabel('反射锥外能量占比 / %'); legend(phaseStateLabels, 'Location', 'best');
title({'雾度评价值（模型代理指标）', ...
    sprintf('圆锥半角 %.1f° / 全角 %.1f°', cfg.metricAngleDeg, 2*cfg.metricAngleDeg)});
grid on; ylim([0, 100]);
angularFigure = figure('Name', 'COE 累计角能量与观察角域', 'Color', 'w', 'Position', [140, 140, 1150, 470]);
angularLayout = tiledlayout(1, 2, 'TileSpacing', 'compact');
title(angularLayout, sprintf('累计绝对反射率 · 中心 (%.2f°, %.2f°)', observation.centerDeg));
nexttile;
plot(withoutAngular.radiusDeg, 100*withoutAngular.radiantCumulativeReflectance, '-', 'LineWidth', 1.5);
hold on;
plot(withAngular.radiusDeg, 100*withAngular.radiantCumulativeReflectance, '--', 'LineWidth', 1.5);
formatCumulativeAxes(observation, phaseStateLabels, spectrumLabel, displayStyle);
nexttile;
plot(withoutAngular.radiusDeg, 100*withoutAngular.photopicCumulativeReflectance, '-', 'LineWidth', 1.5);
hold on;
plot(withAngular.radiusDeg, 100*withAngular.photopicCumulativeReflectance, '--', 'LineWidth', 1.5);
formatCumulativeAxes(observation, phaseStateLabels, photopicLabel, displayStyle);
% 统一本次输出图的字体、线宽与刻度风格，不修改 MATLAB 全局图形默认值。
for outputFigure = [geometryFigure, diffractionFigure, metricsFigure, angularFigure, colorFigure]
    styleCoeFigure(outputFigure, displayStyle);
end

%% 7. 保存物理数据、参数、指标表和图片
% 接收面长度坐标图复用角度图的共同归一化 dB 数据，方便直接对照。
receiverFigure = figure('Name', 'COE 接收面辐照度', 'Color', 'w', ...
    'Position', [100, 80, 1200, 650]);
receiverLayout = tiledlayout(2, numel(previewIndices), 'TileSpacing', 'compact');
title(receiverLayout, sprintf('接收面 z=%.4g m；各波长独立入射 %.4g W/m²', ...
    cfg.propagationDistance, cfg.incidentIrradiance));
subtitle(receiverLayout, sprintf('与角度图共用参考和 dB 色标；I_{ref}=%.4g W/m²', ...
    displayData.referenceIrradiance));
for state = 1:2
    for column = 1:numel(previewIndices)
        k = previewIndices(column);
        nexttile;
        imagesc(result.receiverCoordinate*1e3, result.receiverCoordinate*1e3, ...
            stateDbImages{state}(:, :, column));
        axis image xy; colormap(gca, gray(256)); clim([displayFloorDb, 0]);
        cb = colorbar; cb.Ticks = linspace(displayFloorDb, 0, 4); cb.Label.String = '相对强度 / dB';
        xlabel('X / mm'); ylabel('Y / mm');
        title(sprintf('%s · %.0f nm', phaseStateLabels(state), cfg.wavelengths(k)*1e9));
    end
end
styleCoeFigure(receiverFigure, displayStyle);
if saveResults
    if ~isfolder(outputDir)
        mkdir(outputDir);
    end
    if useProjectConfig
        writelines(projectConfig.text, fullfile(outputDir, caseName + "_config.txt"), 'Encoding', 'UTF-8');
    end
    % MAT 保留未经显示归一化的数组；cfg 记录本次实际计算参数。
    save(fullfile(outputDir, caseName + "_results.mat"), 'cfg', 'result', 'baseline', ...
        'projectName', 'projectConfig', 'inputDir', 'inputFile', 'displayFloorDb', 'displayData', 'comparisonMetrics', 'phaseStateLabels', ...
        'displayStyle', 'useHeightExample', 'heightExample', 'spectrumMode', 'spectrum', ...
        'observation', 'withoutAngular', 'withAngular', 'observationMetrics', 'cumulativeMetrics', ...
        'colorData', 'colorOptions', 'colorView', 'spotOptions', 'diffractionSpots', 'spotMetrics', ...
        'opticalModel', 'stack', 'stackGeometry', '-v7.3');
    if opticalModel == "multilayer"
        writetable(cfg.coeStack.regionTable, fullfile(outputDir, caseName + "_stack_regions.csv"));
        writetable(cfg.coeStack.materialTable, fullfile(outputDir, caseName + "_stack_materials.csv"));
        anodeMask = stackGeometry.anodeMask; bmOpenMask = stackGeometry.bmOpenMask;
        cfRMask = stackGeometry.cfRMask; cfGMask = stackGeometry.cfGMask; cfBMask = stackGeometry.cfBMask;
        save(fullfile(outputDir, caseName + "_layer_patterns.mat"), ...
            'anodeMask', 'bmOpenMask', 'cfRMask', 'cfGMask', 'cfBMask');
    end
    writetable(result.metrics, fullfile(outputDir, caseName + "_metrics.csv"));
    writetable(baseline.metrics, fullfile(outputDir, caseName + "_baseline_metrics.csv"));
    writetable(comparisonMetrics, fullfile(outputDir, caseName + "_comparison_metrics.csv"), 'Encoding', 'UTF-8');
    writetable(observationMetrics, fullfile(outputDir, caseName + "_observation_metrics.csv"), 'Encoding', 'UTF-8');
    writetable(cumulativeMetrics, fullfile(outputDir, caseName + "_cumulative_energy.csv"));
    writetable(spotMetrics, fullfile(outputDir, caseName + "_diffraction_spots.csv"), 'Encoding', 'UTF-8');
    if spectrumMode == "D65"
        exportgraphics(colorFigure, fullfile(outputDir, caseName + "_d65_color.png"), 'Resolution', 180);
        imwrite(colorData.withoutHeight.sRGB, fullfile(outputDir, caseName + "_srgb_linear_without_height_phase.png"));
        imwrite(colorData.withHeight.sRGB, fullfile(outputDir, caseName + "_srgb_linear_with_height_phase.png"));
        spectrumTable = table(spectrum.wavelengthNm, spectrum.d65, spectrum.quadratureNm, ...
            spectrum.cmfXYZ(:, 1), spectrum.cmfXYZ(:, 2), spectrum.cmfXYZ(:, 3), ...
            'VariableNames', {'WavelengthNm', 'D65SPD', 'QuadratureNm', 'XBar', 'YBar', 'ZBar'});
        writetable(spectrumTable, fullfile(outputDir, caseName + "_spectrum.csv"));
    end
    if hasRGB
        imwrite(displayData.withoutHeightRGB, fullfile(outputDir, caseName + "_rgb_without_height_phase.png"));
        imwrite(displayData.withHeightRGB, fullfile(outputDir, caseName + "_rgb_with_height_phase.png"));
    end
    exportgraphics(geometryFigure, fullfile(outputDir, caseName + "_geometry.png"), 'Resolution', 180);
    exportgraphics(diffractionFigure, fullfile(outputDir, caseName + "_diffraction.png"), 'Resolution', 180);
    exportgraphics(receiverFigure, fullfile(outputDir, caseName + "_receiver_irradiance.png"), 'Resolution', 180);
    exportgraphics(metricsFigure, fullfile(outputDir, caseName + "_metrics.png"), 'Resolution', 180);
    exportgraphics(angularFigure, fullfile(outputDir, caseName + "_cumulative_energy.png"), 'Resolution', 180);
    fprintf('计算完成，结果已保存到：%s\n', outputDir);
end
% 运行结束后 cfg、result、baseline 保留在工作区，方便继续查看与分析。
