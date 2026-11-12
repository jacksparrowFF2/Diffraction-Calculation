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
inputDir = fullfile(moduleDir, 'Input');
outputDir = fullfile(moduleDir, 'Output');
inputFile = "";                  % Input/ 下的 MAT 文件名，如 "design.mat"；空字符串表示生成默认结构
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
cfg.radius = 10e-6;             % 椭圆短半轴；aspectRatio=1 时为圆半径
cfg.aspectRatio = 1;            % 长半轴 / 短半轴
cfg.cellAnglesDeg = [30, 0; 0, -30]; % 2×2 晶胞内各单元的旋转角；圆形时旋转不改变形状
cfg.periods = [16, 16];         % [x 方向周期数, y 方向周期数]，均为偶数
cfg.samplesPerPeriod = 32;      % 每周期采样点数；采样间距 = pitch / samplesPerPeriod
cfg.paddingFactor = 4;          % FFT 补零倍数，影响远场采样密度
cfg.pupilRadius = Inf;          % Inf 表示整个矩形计算面参与反射
cfg.mask = [];                 % 留空时自动生成几何；自定义矩阵尺寸见下方 gridSize

% 2.2 光谱与反射面
spectrumMode = "D65";          % "D65"：380~780 nm 色度合成；"RGB"：原三波长演示
wavelengthStepNm = 5;          % 光谱步长；可改为 1/2/5/10 nm，越细计算越久
previewWavelengthNm = [620, 550, 450]; % 远场图仅预览这些波长，全部波长均参与积分
cfg.wavelengths = [620, 550, 450]*1e-9; % RGB 模式下的三个代理波长
cfg.spectralWeights = [1, 1, 1]; % RGB 模式使用的入射权重
cfg.channelResponse = eye(3);   % RGB 模式的 R/G/B 通道
mustBeMember(spectrumMode, ["D65", "RGB"]);
if spectrumMode == "D65"
    validateattributes(wavelengthStepNm, {'numeric'}, {'scalar', 'integer', 'positive', '<=', 20});
    assert(mod(400, wavelengthStepNm) == 0, 'COEReflection:SpectrumStep', '波长步长必须整除 400 nm。');
    spectrum = loadCOESpectrum(380:wavelengthStepNm:780);
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
stack.geometry.source = "diamond"; % diamond：可调示例；mat：同一 MAT；images：独立图片
stack.geometry.matFile = "coe_layers.mat"; % 字段 anodeMask、cfR/G/BMask，可选 bmOpenMask 或 bmMask
stack.geometry.anodeFile = "anode.png";
stack.geometry.bmFile = "";     % 空时自动从阳极外扩生成；白色默认代表开口
stack.geometry.bmImageIsOpening = true; % false 表示 BM 图片白色是实心材料
stack.geometry.cfFiles = ["CF_R.png", "CF_G.png", "CF_B.png"];
stack.bmExpansion = 5e-6;       % BM 开口相对阳极边界向外扩张 5 μm，可调
stack.demo.anodeHalfDiagonal = [12, 10, 12]*1e-6; % R/G/B 菱形半对角线，可调
stack.demo.cfExpansion = 7e-6;  % CF 比 BM 开口更大；各色不得交叠
stack.demo.colorTile = [1, 2; 2, 3]; % R/G/G/B 演示排布；Pitch 由 cfg.pitch 定义为相邻单元间距
stack.anodeReflectance = 1.0;   % 仅阳极图样内的裸底为 100% 反射
stack.anodePhaseDeg = 180;      % 理想电反射边界相位；阳极外为匹配吸收底（裸底 r=0）
stack.materials.BM.n = 1.8; stack.materials.BM.k = 0.35; stack.materials.BM.thickness = 1.0e-6;
stack.materials.R.n = 1.65; stack.materials.R.k = 0.03; stack.materials.R.thickness = 1.1e-6;
stack.materials.G.n = 1.60; stack.materials.G.k = 0.04; stack.materials.G.thickness = 1.0e-6;
stack.materials.B.n = 1.55; stack.materials.B.k = 0.05; stack.materials.B.thickness = 0.9e-6;
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
    fprintf('使用脚本参数生成默认反射面。\n');
    if useHeightExample
        cfg.heightMap = createCOEHeightExample(cfg, heightExample); %#ok<UNRCH> 可选示例默认关闭
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
    [cfg.coeStack, stackGeometry] = prepareCOEStack(cfg, stack, inputDir);
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
withoutAngular = computeCOEAngularMetrics(baseline, observation, spectrum);
withAngular = computeCOEAngularMetrics(result, observation, spectrum);
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
diffractionSpots = computeCOEDiffractionSpots(baseline, result, spectrum, spotOptions);
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
        synthesizeCOEColor(baseline, result, spectrum, colorOptions);
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
        composeCOEReflectionRGB(baseline.channelSignal, result.channelSignal, rgbOptions);
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
columnCount = numel(previewIndices) + double(hasRGB);
diffractionLayout = tiledlayout(2, columnCount, ...
    'TileSpacing', 'compact', 'Padding', 'compact');
title(diffractionLayout, '同一结构 · 高度相位作用对比');
annotationHint = "辅助标记已关闭";
if displayStyle.showAnnotations
    annotationHint = "辅助标记按 displayStyle 分类开关显示";
end
subtitle(diffractionLayout, spectrumLabel + "｜上/下：两种高度相位状态；" + annotationHint);
stateResults = {baseline, result};
stateImages = {};
if hasRGB
    stateImages = {displayData.withoutHeightRGB, displayData.withHeightRGB};
end
for state = 1:2
    for column = 1:numel(previewIndices)
        k = previewIndices(column);
        nexttile((state-1)*columnCount + column);
        imagesc(angleDeg, angleDeg, 10*log10(max( ...
            stateResults{state}.spectralEnergy(:, :, k)/referenceEnergy, 10^(displayFloorDb/10))));
        formatFarFieldAxes(gca, cfg, displayStyle);
        colormap(gca, gray(256)); clim([displayFloorDb, 0]);
        cb = colorbar; cb.Ticks = linspace(displayFloorDb, 0, 4); cb.Label.String = '相对能量 / dB';
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
    styleCOEFigure(outputFigure, displayStyle);
end

%% 7. 保存物理数据、参数、指标表和图片
if saveResults
    if ~isfolder(outputDir)
        mkdir(outputDir);
    end
    % MAT 保留未经显示归一化的数组；cfg 记录本次实际计算参数。
    save(fullfile(outputDir, caseName + "_results.mat"), 'cfg', 'result', 'baseline', ...
        'inputFile', 'displayFloorDb', 'displayData', 'comparisonMetrics', 'phaseStateLabels', ...
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
    exportgraphics(metricsFigure, fullfile(outputDir, caseName + "_metrics.png"), 'Resolution', 180);
    exportgraphics(angularFigure, fullfile(outputDir, caseName + "_cumulative_energy.png"), 'Resolution', 180);
    fprintf('计算完成，结果已保存到：%s\n', outputDir);
end
% 运行结束后 cfg、result、baseline 保留在工作区，方便继续查看与分析。

function formatFarFieldAxes(ax, cfg, style)
%formatFarFieldAxes 统一角度坐标，并标出雾度评价圆锥边界。
    axis(ax, 'image');
    ax.YDir = 'normal';
    limit = cfg.angleHalfRangeDeg;
    xlim(ax, [-limit, limit]); ylim(ax, [-limit, limit]);
    ax.XTick = linspace(-limit, limit, 5);
    ax.YTick = ax.XTick;
    xlabel(ax, 'θ_x / °'); ylabel(ax, 'θ_y / °');
    if ~style.showAnnotations
        return;
    end
    hold(ax, 'on');
    phi = linspace(0, 2*pi, 361);
    if style.showMetricCone
        plot(ax, cfg.metricAngleDeg*cos(phi), cfg.metricAngleDeg*sin(phi), ...
            '--', 'Color', [0.6, 0.85, 0.9], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
    obs = cfg.observation;
    if style.showObservationRegion
        plot(ax, obs.centerDeg(1) + obs.halfAngleDeg*cos(phi), obs.centerDeg(2) + obs.halfAngleDeg*sin(phi), ...
            ':', 'Color', [1, 0.7, 0.2], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
    if style.showObservationRegion && obs.innerHalfAngleDeg > 0
        plot(ax, obs.centerDeg(1) + obs.innerHalfAngleDeg*cos(phi), ...
            obs.centerDeg(2) + obs.innerHalfAngleDeg*sin(phi), ...
            ':', 'Color', [1, 0.7, 0.2], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
end

function formatCumulativeAxes(observation, labels, panelTitle, style)
%formatCumulativeAxes 半径内累计值始终以入射量为分母，不强制将曲线末端归一化到 100%。
    if style.showAnnotations && style.showObservationRegion
        xline(observation.halfAngleDeg, ':', '观察外半角', 'HandleVisibility', 'off');
    end
    if style.showAnnotations && style.showObservationRegion && observation.innerHalfAngleDeg > 0
        xline(observation.innerHalfAngleDeg, ':', '观察内半角', 'HandleVisibility', 'off');
    end
    xlabel('相对观察中心的累计半角 / °'); ylabel('累计绝对反射率 / %');
    title(panelTitle); legend(labels, 'Location', 'best'); grid on;
    xlim([0, max(xlim)]);
end

function markDiffractionSpots(ax, spots, groupName, wavelength, style)
%markDiffractionSpots 标出共同积分圆、编号与零级排除圆；不使用显示增强结果找峰。
    if ~style.showAnnotations
        return;
    end
    groupIndex = find(string({spots.groups.name}) == groupName);
    if groupName == "wavelength"
        groupIndex = groupIndex(abs([spots.groups(groupIndex).wavelengthNm] - wavelength) < 1e-6);
    end
    if isempty(groupIndex)
        return;
    end
    centers = spots.groups(groupIndex(1)).centersDeg;
    phi = linspace(0, 2*pi, 181);
    hold(ax, 'on');
    if style.showZeroOrderExclusion
        plot(ax, spots.options.zeroOrderExclusionDeg*cos(phi), ...
            spots.options.zeroOrderExclusionDeg*sin(phi), ':', 'Color', [1, 0.35, 0.35], ...
            'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
    for j = 1:size(centers, 1)
        if style.showSpotRegions
            plot(ax, centers(j, 1) + spots.options.integrationHalfAngleDeg*cos(phi), ...
                centers(j, 2) + spots.options.integrationHalfAngleDeg*sin(phi), ...
                '-', 'Color', [1, 0.35, 0.8], 'LineWidth', 0.8, 'HandleVisibility', 'off');
        end
        if ~style.showSpotLabels
            continue;
        end
        limitsX = xlim(ax); limitsY = ylim(ax);
        if centers(j, 1) >= limitsX(1) && centers(j, 1) <= limitsX(2) ...
                && centers(j, 2) >= limitsY(1) && centers(j, 2) <= limitsY(2)
            direction = centers(j, :)/max(norm(centers(j, :)), eps);
            span = max(diff(limitsX), diff(limitsY));
            offset = max(1.3*spots.options.integrationHalfAngleDeg, 0.065*span);
            labelPosition = centers(j, :) + offset*direction;
            labelPosition(1) = min(max(labelPosition(1), limitsX(1) + 0.04*span), limitsX(2) - 0.04*span);
            labelPosition(2) = min(max(labelPosition(2), limitsY(1) + 0.04*span), limitsY(2) - 0.04*span);
            if offset > 1.5*spots.options.integrationHalfAngleDeg
                start = centers(j, :) + spots.options.integrationHalfAngleDeg*direction;
                plot(ax, [start(1), labelPosition(1)], [start(2), labelPosition(2)], ...
                    '-', 'Color', [0.8, 0.3, 0.65], 'LineWidth', 0.5, 'HandleVisibility', 'off');
            end
            text(ax, labelPosition(1), labelPosition(2), "S" + j, ...
                'Color', [1, 0.75, 0.9], 'HorizontalAlignment', 'center', ...
                'FontSize', 8.5, 'Tag', 'COESpotLabel', 'Clipping', 'on');
        end
    end
end

function styleCOEFigure(fig, style)
%styleCOEFigure 统一中文字体、坐标轴和色标样式，仅作用于本次生成的图。
    set(findall(fig, '-property', 'FontName'), 'FontName', style.fontName);
    set(findall(fig, '-property', 'FontSize'), 'FontSize', style.fontSize);
    set(findall(fig, 'Type', 'axes'), 'LineWidth', style.lineWidth, ...
        'TickDir', 'out', 'Box', 'on', 'TickLabelInterpreter', 'none');
    set(findall(fig, 'Type', 'colorbar'), 'LineWidth', style.lineWidth, 'TickDirection', 'out');
    set(findall(fig, 'Type', 'legend'), 'Box', 'off');
    set(findall(fig, 'Tag', 'COESpotLabel'), 'FontSize', 8.5);
end
