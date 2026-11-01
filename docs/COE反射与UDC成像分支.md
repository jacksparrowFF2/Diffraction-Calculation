# COE 反射细化与 UDC 成像分支

开发分支：`feature/coe-reflection-udc`。原有脚本及参考代码均保留。

## 直接运行

在项目根目录的 MATLAB 命令窗口：

```matlab
addpath('src');
reflection = runCOEReflection();
[udc, images] = runUDCImaging();
demoCOEUDC();
```

演示将概览 PNG、反射/透射指标 CSV 和包含完整参数、PSF、线性图像的 MAT 文件写入 `ImageForShow/COE_UDC/`。该目录沿用原项目的忽略规则，不进入 Git。运行不依赖当前文件夹中的输入图片，默认测试图路径由函数自身位置定位。建议 MATLAB R2020b 以上；本次验证使用 R2024a。新增计算只用基础 MATLAB，不需要图像处理工具箱。

## 与原程序的关系

- 沿用 `CustomAperture_fraunhofer_fft.m` 的 50 µm 周期、10 µm 半径、2×2 晶胞和 `CustomPattern` 形状函数。`aspectRatio` 非 1 时才能看出四孔旋转差异。
- 沿用 `CustomApertureWithPhase_fraunhofer_fft.m` 中反射双程高度相位的定义。现在几何、高度图与材料系数分离，同一高度图也可以传给透射分支。
- 图像卷积复用已有 `fftconv2_same.m`，使用奇数 PSF 核；增加边缘延拓，避免黑色填充造成的边缘变暗。
- 新核心使用明确物理采样的 FFT，再把**能量**按面积重叠积分到统一像面。原 `fraunhofer_fft_SPH` / `myFFT2_SPH` 仍保留，便于复查原结果。
- 原参考代码 `WangAnPing_code/衍射/*spectrum-*.m` 将 D65 强度乘反射强度后赋给 `U`，再取 `abs(FFT(U))^2`，会再次平方这些强度权重；原主脚本图像模式只用了末个波长的 `I`，并按峰值归一化。新增入口通过复振幅建场、逐波长强度相加和能量归一化处理这两点。

## 物理模型与适用范围

这是法向入射、标量、薄元件、近轴模型。反射计算远场角分布；UDC 假设屏幕与镜头入瞳共面、远景对焦、视场内空间不变。屏幕至入瞳的非零间距、有限物距、离轴入射、有限扩展光源、多层界面/腔体反射和偏振均未建模。UDC 不是把人眼观察距离改成镜头焦距后继续做自由空间传播，而是使用入瞳复透射场的焦面傅里叶变换。

默认的反射率、透过率、入瞳尺寸、相位高度都是示例参数，不代表具体屏幕。R/G/B 采用 620/550/450 nm 三个窄带代理；不是 D65 色度积分，也没有相机实测光谱响应。MTF 是有限窗口、像素积分后的系统代理，不等于相机实测分辨率或画质评价。

反射复场：

```text
E_r(x,y,λ) = P(x,y) · r(x,y,λ) · exp[i · 4π · Δn · h(x,y)/λ]
```

透射复场：

```text
E_t(x,y,λ) = P(x,y) · t(x,y,λ) · exp[i · 2π · Δn · h(x,y)/λ]
```

强度系数 R/T 使用平方根转换为振幅；复数 `r/t` 可以直接输入。各波长分别计算 `|FFT(E)|²` 后非相干相加。高度相位是法向、简化光程差模型；如复数材料系数已包含该高度相位，应把 `heightMap=0`，避免重复计入。

## 参数与自定义几何

```matlab
cfg = displayDiffractionDefaults("reflection");
cfg.featureIntensity = [0.55 0.60 0.62]; % 强度反射率 R，依次对应 wavelengths
cfg.backgroundIntensity = 0.01;
cfg.aspectRatio = 1.3;
cfg.heightMap = 220e-9; % 均匀高度只产生整体相位，不会改变衍射强度
reflection = runCOEReflection(cfg);
```

默认 `periods=[16 16]` 是 `[Nx Ny]`，每周期 32 个采样，得到 512×512 的物面；周期端点不重复。`paddingFactor=4` 得到 2048×2048 FFT。非方形周期数也支持。两方向周期数必须为偶数。

如已有旧脚本生成的 `BaseMask`、`h_map`，按其采样与物理尺寸设置：

```matlab
cfg = displayDiffractionDefaults("reflection");
cfg.pitch = Pitch;
cfg.periods = [Nx_period Ny_period];
cfg.samplesPerPeriod = nn;
cfg.mask = BaseMask;
cfg.heightMap = h_map;
cfg.deltaIndex = n_high - n_low;
cfg.paddingFactor = 1; % 大网格先控制内存，再检查收敛
reflection = runCOEReflection(cfg);
```

`mask` 和 `heightMap` 必须与 `[Ny*nn, Nx*nn]` 网格一致。旧脚本用 `linspace` 含两端点；严格定量比较时，应在新网格重新生成或按物理坐标重采样几何/高度，不能仅因数组同尺寸就认为坐标完全一致。特征小于采样间距时应增加 `samplesPerPeriod`，补零不能恢复丢失的结构。

默认 `mask=1` 的区域代表原几何孔的特征区域，`featureIntensity/backgroundIntensity` 由调用者指定其物理含义；反射区域和透光区域不必相同。真实器件建议分别给出两分支的复数图。

```matlab
cfg.complexCoefficient = rMap; % reflection: 复反射系数，支持二维或 Ny×Nx×L
cfg.heightMap = 0;             % 当 rMap 已含完整相位时
```

`complexCoefficient` 非空时取代 feature/background 系数，但仍乘入瞳和额外高度相位。要求被动、功率归一化后的 `abs(coefficient)<=1`。透射到不同折射率介质时，需要先将电场 t 按介质导纳/功率通量关系转换为功率归一化系数，不能直接把电场 t 的模平方当透过率。复振幅图通常优于从显微照片直接当作物理透过率。

## 反射输出

```matlab
cfg = displayDiffractionDefaults("reflection");
cfg.angleHalfRangeDeg = 6; % 近轴角坐标 ±6°，限制到 ±10°以内
cfg.metricAngleDeg = 2.5;
r = runCOEReflection(cfg);
disp(r.metrics);
```

`coordinate` 为近轴角度弧度，满足 `θ≈λ·空间频率`，并非高角度精确球面坐标。主要指标：

| 字段 | 定义 |
|---|---|
| Throughput | 全部离散傅里叶平面能量 / 理想入射入瞳能量，即平均 R 或 T |
| CapturedThroughput | 输出窗口中的能量 / 理想入射能量 |
| CapturedFraction | 输出窗口能量 / 全部计算能量 |
| ParsevalAbsoluteError | FFT 能量与物面能量之差的绝对值 |
| ZeroOrderFraction | 半径 `0.45·λ/Pitch` 内的能量 / 全部反射能量；这是中心区域判据 |
| OutsideAngleFraction | `1 - 指定角内能量/全部计算能量`，含输出窗口外尾部 |
| ZeroOrderRetention | 有高度图与无高度图的零级区域能量比例；分母为零时为 NaN |

`baselineKernel/baselineMetrics` 是清除额外高度图后的同材料、同几何基准；`complexCoefficient` 内已有的相位仍保留。`OutsideAngleFraction` 仅为模型角域能量指标，不能称为标准仪器雾度。所有指标依赖物面、FFT 和输出角度网格收敛；FFT 高频部分仍采用标量近轴延拓，不解释为真实高角度/倏逝波能量。

## UDC 图像与光谱

```matlab
cfg = displayDiffractionDefaults("udc");
cfg.focalLength = 3e-3;
cfg.pupilRadius = 0.35e-3;
cfg.sensorPitch = 1e-6;
cfg.outputSize = 401;
[u, im] = runUDCImaging(cfg); % 默认读取项目 USAF500.png 并解码 sRGB

rgb = imread('my_scene.png');
[u, im] = runUDCImaging(cfg, rgb, "srgb");
hdr = zeros(501, 501, 3);
hdr(251, 251, :) = 1000;
im = formUDCImage(hdr, u, "linear"); % 不提前裁剪 HDR，也不重新计算光学 PSF
```

输入解释为**已映射到传感器面的理想场景辐照度**，其采样间距为 `sensorPitch`。程序不从照片像素自动推导物距、FOV 或倍率；普通照片仅作为示意测试目标。默认 0.7 mm 入瞳与 3 mm 焦距约为 f/4.3，是满足当前几何窗口的演示设置。若修改入瞳，物面宽高都必须覆盖完整直径；不能让数组边界静默切掉镜头孔径。

- `kernel`：各通道实际能量核，和为 `channelCapturedThroughput`，含材料损耗和有限核损失。
- `psf`：各非零核按**总和**归一到 1，表示模糊形状；零透射通道保持全零。
- `clearKernel/clearPsf`：同一入瞳去掉屏幕后的镜头衍射基准。
- `images.clear`：无屏幕、同曝光的有限核结果。
- `images.fixedExposure`：有屏幕、同曝光结果，不消除变暗。
- `images.brightnessMatched`：单位和核结果，只比较模糊。不能解释成数字增益后的噪声效果，程序目前没有传感器噪声模型。
- `mtf/spatialFrequency`：每通道二维归一化 MTF 和 cycles/mm 频率轴。无信号通道 MTF 为 0。UDC 不使用反射角指标，对应表列为 NaN。

卷积在线性强度域进行、各通道独立、默认边缘重复延拓。保存的线性结果保留 HDR；仅演示 PNG 编码时裁剪到显示范围。原始场景视野之外的亮点不会凭空生成，需要扩大输入画布模拟。

支持更密波长网格。更改波长时，必须同步提供权重和波长×通道响应：

```matlab
cfg.wavelengths = (450:5:650)*1e-9;
cfg.spectralWeights = sourceSPD(:).'*5; % 例：每 nm 的光谱 × 5 nm 积分宽度
cfg.channelResponse = measuredResponse; % L×3 非负相机响应，行与波长一一对应
```

每通道 `kernel` 除以本通道入射光谱积分，所以是相对该通道无损入射的能量比；`channelSignal` 保留权重未归一化的积分，用于观察光源谱/响应造成的通道相对强弱。默认不含 CIE CMF、色适应或相机颜色矩阵，不能把其三通道直接视为预测实物颜色。

## 数值检查与收敛

有限 PSF 窗口截获不足 98% 时会警告；原生 FFT 像面间距大于传感器像元时也会警告。前者通过增大 `outputSize` 检查，后者通过增加 `paddingFactor` 检查。物面采样通过增加 `samplesPerPeriod` 检查；几何尺寸不应随采样变化。窄衍射峰、亚像元结构、高度图边界和细线金属需要特别检查。重采样按每个 FFT 栅格内分段常量能量密度处理；补零和网格加密收敛后才适合定量结论。

## 自动测试

在项目根目录运行：

```matlab
results = runtests('tests/DisplayDiffractionTest.m');
assertSuccess(results);
```

测试类自动加入并恢复 `src` 路径，使用基础 MATLAB 的 `matlab.unittest`。21 项测试覆盖已知反射/透射强度、Parseval 能量守恒、圆形入瞳面积加权、同一高度图单程 π/2 与双程 π 的解析零级结果、均匀高度整体相位不变性、正一阶角度/传感器位置、光谱非相干强度叠加、RGB 独立卷积与空间域 `conv2` 对照、HDR 脉冲幅度与方向、sRGB 解码、零反射/透射边界、有限窗口损失与警告、MTF 归一化、入瞳截断及有源系数输入拒绝。

解析测试使用小型矩形入瞳，输出格与 FFT 格对齐，分别验证 `R=0.36`、`T=0.25`，以及相位光栅的一阶位置 `λ/Pitch` / `f·λ/Pitch`。这些测试验证计算一致性，不代替真实器件或完整连续模型的收敛和实测验证。默认完整示例另行在 R2024a 运行并检查预览；测试结果汇总保存在 `ImageForShow/COE_UDC/test_results.csv`。

测试也修正了相等采样间距受浮点舍入触发“采样不足”误报的问题：比较时使用 `1e-12` 相对容差，并设专门的无警告回归用例。

## 模型依据

- [MathWorks：二维傅里叶变换与孔径衍射](https://www.mathworks.com/help/matlab/math/two-dimensional-fft.html)
- [Zhou 等，Image Restoration for Under-Display Camera，CVPR 2021](https://openaccess.thecvf.com/content/CVPR2021/html/Zhou_Image_Restoration_for_Under-Display_Camera_CVPR_2021_paper.html)：屏幕结构、PSF 与 UDC 图像形成的研究背景。
- [Feng 等，Removing Diffraction Image Artifacts in Under-Display Camera，CVPR 2021](https://openaccess.thecvf.com/content/CVPR2021/html/Feng_Removing_Diffraction_Image_Artifacts_in_Under-Display_Camera_via_Dynamic_Skip_CVPR_2021_paper.html)：强光衍射伪影和物理成像建模背景。这里未实现论文的恢复网络。
