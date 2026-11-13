# COE 反射衍射

当前默认从 `Input/项目名称/config.txt` 读取参数，配置值优先于主脚本默认值。正式使用只需选择项目并编辑其配置；不调整像素形状。语法见 [项目配置说明](项目配置说明.md)。

实际版图和 nk 数据的准备、主脚本修改示例见 [实际图样与材料输入指南](实际图样与材料输入指南.md)。

本模块是后续主要开发方向之一，用于评估产品 COE 反射设计方案的可行性，并在一致的计算条件与评价指标下对比不同方案的效果。

本目录可独立复制，`src/private/` 随目录一起保留。

操作和迁移请优先阅读 [程序使用说明](程序使用说明.md)；继续开发请阅读 [项目交接文档](项目交接文档.md)。另一台电脑直接复制完整模块文件夹或从 Git 拉取，无需迁移 ZIP。

## 主脚本：从输入到结果保存

在 MATLAB 中打开并运行 `src/COEReflection_main.m`。脚本按以下顺序展开，可从上往下阅读和修改：

0. **初始化环境**：执行 `clearvars`、`clc`、`close all`，保留 MATLAB 路径和 MCP 连接。
1. **文件输入与输出设置**：指定 `Input/` 下的 MAT 文件名、输出前缀和保存开关。
2. **模型参数**：集中设置几何、采样、光谱、反射率、高度和观察角度，附中文说明及单位。
3. **读取设计文件并检查尺寸**：将输入的反射面数据应用到本次计算参数。
4. **叠加高度相位**：调用底层函数计算双程相位、FFT 和远场能量。
5. **不叠加高度相位及指标对比**：保持其它参数不变，仅清除额外高度图，计算零级保持率。
6. **结果展示**：材料区域和高度图、选定波长及 D65 色度合成图、雾度评价图、累计角能量曲线。
7. **结果保存**：输出实际参数、物理数组、CSV 指标表和 PNG 图片。

默认 `opticalModel="multilayer"`，无需准备文件即可运行可调 Diamond 示例；采用 D65 380–780 nm、5 nm 步长（81 个波长），逐位置计算阳极/BM/RGB CF 的相干膜堆复反射系数。默认 `cfg.heightMap = 0`（0 nm）、`useHeightExample = false`，暂不添加额外高度相位，两种状态的结果应一致。膜厚和不同顶面到共同参考面的传播相位始终保留。每次运行会清理工作区和旧图窗。

### 阳极 / BM / RGB CF 多层输入

主脚本的 `stack` 参数集中定义图样来源、阳极反射边界、BM 外扩、Diamond 尺寸，以及 BM/R/G/B 的 n、k、厚度。详细定义见 [简化 COE 多层模型约定](简化COE多层模型约定.md)。

- 正式仿真只导入项目图片或 MAT 图样，不提供像素尺寸、圆角、倾角等形状调整。预生成示例按参考 PDL 轮廓构建阳极，R/B 交错、G 位于四邻中心且交替倾斜，R:G:B 数量比 1:2:1。
- BM 开口默认由阳极边界按 5 μm 的离散 Euclidean 距离外扩，CF 按 7 μm 外扩并校验不交叠。阳极外环带是吸收底，绝不当作扩大的阳极。
- `source="mat"`：读取 `Input/coe_layers.mat`，字段 `anodeMask`、`cfRMask`、`cfGMask`、`cfBMask`，可选 `bmOpenMask`（1=开口）或 `bmMask`（1=材料）。不提供 BM 时自动外扩。
- `source="images"`：分别指定阳极、BM、R/G/B CF 的图片文件。白色表示对应区域，BM 的白色含义由 `bmImageIsOpening` 指定；BM 文件名为空时自动外扩。图像不自动缩放，所有图样须与计算网格同尺寸且对齐。
- `stack.materials.BM/R/G/B` 各有 `n`、`k`、`thickness`（m）。n/k 支持常数或逐波长向量；非空 `nkFile` 优先读取三列 CSV `[波长 nm,n,k]`，必须覆盖整个计算波段。
- 阳极裸底默认 `R=1`、相位 180°。阳极外默认裸底 `r=0`，按与入射介质阻抗匹配的吸收终端处理；覆盖膜层后仍计算界面反射和干涉，最终局部反射不一定为零。
- 输出增加 `_stack_regions.csv`、`_stack_materials.csv` 与 `_layer_patterns.mat`。最后一个文件可以复制到 `Input/`，切到 MAT 模式复用。

若需复现旧模型，设置 `opticalModel="simple"`。旧参数 `cfg.radius/aspectRatio/cellAnglesDeg/materials` 和旧单一 `mask` 仅在该模式下使用。外部 `complexCoefficient` 非空时仍具有最高优先级。

### 光谱与材质

- `spectrumMode = "D65"`：读取本目录 `Input/Colorimetry/` 的 CIE 官方 D65 和 CIE 1931 2°标准观察者数据。`wavelengthStepNm` 可设为 1/2/5/10 nm 等整除 400 的步长。默认 5 nm 是计算成本与精度的初始折中，最终方案应检查波长和空间采样收敛。
- `previewWavelengthNm` 只控制远场对比图显示哪些单波长行，不减少参与光谱积分的波长数。
- 旧 `simple` 模式中，`cfg.materials.enabled = true` 时，mask=1/0 分别代表两种单界面材料；其常数示例不对应指定产品。
- 旧单界面系数为 `r=(n_incident-(n+i*k))/(n_incident+(n+i*k))`，`R=abs(r)^2`；多层模式使用独立膜堆递推。额外高度相位仍使用单独的有效 `deltaIndex`，不将材料 n 自动当成该光程差。
- 输入 `complexCoefficient` 非空时优先于材料 n/k；关闭材料模型时回到 `featureIntensity/backgroundIntensity` 的强度反射率模型。
- `spectrumMode = "RGB"` 保留原三波长增强预览。复现早先 0.60/0.01 强度模型，还需设置 `opticalModel="simple"`、`cfg.materials.enabled=false`。该模式的第二条累计曲线为代理权重，不代表标准光度结果。

### 观察角域与累计角能量

在主脚本中设置：

```matlab
observation.centerDeg = [0, 0];       % [theta_x, theta_y] 中心方向
observation.halfAngleDeg = 2.5;       % 观察外半角
observation.innerHalfAngleDeg = 0;    % 0=圆锥；正值=排除中心后的环带
observation.curveSamples = 301;
```

角域必须完整落在输出窗口内，程序会拒绝越界区域，避免静默少算能量。该设置独立于 `cfg.metricAngleDeg` 的雾度评价阈值。

观察角域绝对反射率为 `E_view/E_incident`；明视觉值为同一区域反射光通量与入射光通量之比。两者分别按 D65 辐射度权重及 D65×y_bar 权重积分，单位 %。

累计曲线表示以观察中心为圆心、半角从 0 逐渐增加的**实心圆锥**内累计绝对反射率；即使观察区域设为环带，累计曲线也不扣除内孔。半径范围只到能完整放入方形输出窗口的最大圆，不强行将末端归一化到 100%。有限角度网格以 bin 中心判断是否落入角域，因此边界附近需检查采样收敛，零半径的累计值定义为 0。

### 非零级衍射斑统计

主脚本的 `spotOptions` 可设置零级排除半角（默认 0.30°）、积分半角（0.15°）、最小中心间距（0.35°）、每组最多光斑数（8）和相对峰阈值（0.08）。`wavelengthNm` 默认统计预览波长，空数组可统计全部波长，另始终计算宽带 `radiant/photopic` 两组。

检测基于两状态原始能量图的共同候选，分别在同一组圆域中积分，输出以入射量为分母的绝对反射率。圆域不进入零级排除范围、不重叠、不越过输出窗口。图中品红圈和 S 编号对应 `_diffraction_spots.csv`。编号不是理论衍射阶次，不同组的能量不能直接相加。

### 主脚本中的非均匀高度示例

高度示例默认关闭。手动设置 `useHeightExample = true` 时，才生成覆盖整个计算面的二值高度条纹：

- 高区为 `heightExample.height = 687.5e-9` m，低区为 0。
- 周期为 `heightExample.period = cfg.pitch/4`，默认 12.5 μm，占空比 50%。
- `heightExample.direction = "x"` 表示沿 x 交替、条纹平行 y；可改为 `"y"`。
- Δn=0.20 时，550 nm 波长的高低区双程相位差为 π；红、蓝分别约 0.887π、1.222π。
- 高度作用于特征区和背景区，但不改变材料区域或强度反射率。此结构仅用于演示，不代表产品或优化结果。

高度周期至少需要 4 个采样点，周期和高区宽度均须为采样间距的整数倍；程序会检查。将 `useHeightExample` 改为 `false` 即恢复手动设置 `cfg.heightMap`（初始值 0）；指定非空 `inputFile` 时跳过示例，使用文件数据及原有参数，不会混入演示高度。

早期三波长、强度反射率配置的量化对比见 [非均匀高度示例与反射评价指标](非均匀高度示例与反射评价指标.md)。当前默认 D65/nk 的说明及结果见 [光谱色度与观察角域](光谱色度与观察角域.md)，两套参数不能直接当作同一次配置比较。

`heightMap` 是空间高度分布 h(x,y)，用于引入空间相位 `4*pi*deltaIndex*h/lambda`，并不是“引入平均相位”。均匀高度仅改变整体相位，不改变远场强度；非均匀高度会重新分配衍射能量。图表统一使用“不叠加高度相位 / 叠加高度相位”，表示同一输入结构的两种计算状态，没有在计算中修改设计方案。`baseline` 是为兼容既有接口保留的变量名。

### 自定义设计输入

MAT 文件支持以下任意字段，文件中的字段覆盖脚本中同名参数，其它参数仍使用脚本设置：

| 字段 | 内容 |
|---|---|
| `mask` | 0~1 材料区域图，推荐二值；1 为特征区，0 为背景区 |
| `heightMap` | 标量或二维高度图，单位 m |
| `complexCoefficient` | 复振幅反射系数 r；标量、二维矩阵或 Ny×Nx×波长数数组，模不大于 1 |

二维矩阵必须覆盖整个计算面，尺寸为 `[periods(2), periods(1)]*samplesPerPeriod`（行 y、列 x），不会自动插值。`mask` 为空时按几何参数自动生成。`complexCoefficient` 非空时优先使用输入系数，替代由 `mask` 和强度反射率构造的系数；额外高度相位仍会叠加。若 r 已包含该高度相位，应将 `heightMap` 设为 0，避免重复计入。“不叠加高度相位”仅清除 `heightMap`，不会清除 r 自身携带的相位，因此不等于“无相位”。

例如，先在 MATLAB 中生成一份与默认采样匹配的非均匀高度输入：

```matlab
heightMap = zeros(1024, 1024); % 默认 16 个周期 × 每周期 64 点
heightMap(:, 257:end) = 200e-9; % 右半面增加 200 nm 高度
save('topic_coe_reflection/Input/design.mat', 'heightMap'); % 在仓库根目录执行
```

然后将主脚本的 `inputFile` 改成 `"design.mat"`，运行整个脚本。

### 输出与显示约定

默认写入 `Output/`，同一 `caseName` 再运行会覆盖同名文件：

- `coe_reflection_results.mat`：`cfg`、`result`（叠加高度相位）、`baseline`（不叠加高度相位）、光谱权重、角域指标、累计曲线、XYZ/xy/线性 sRGB 密度和共同曝光参数。
- `coe_reflection_metrics.csv`、`coe_reflection_baseline_metrics.csv`：分别为叠加 / 不叠加高度相位的逐波长物理指标，比例为 0~1。
- `coe_reflection_comparison_metrics.csv`：含明确状态名、圆锥半角与全角、百分比形式的雾度评价值和零级占比。
- `coe_reflection_observation_metrics.csv`：观察角域的辐射度/明视觉绝对反射率、总反射率、角域中心和内外半角；D65 模式还含该区域积分后的 xy 色坐标。
- `coe_reflection_cumulative_energy.csv`、`coe_reflection_cumulative_energy.png`：两种高度相位状态的累计绝对反射率。
- `coe_reflection_diffraction_spots.csv`：非零级光斑的共同中心、实际峰角位置、积分半角及绝对反射率。
- `coe_reflection_spectrum.csv`：D65 模式下的波长、入射谱、积分步长与三刺激值函数。
- `coe_reflection_geometry.png`、`coe_reflection_diffraction.png`、`coe_reflection_metrics.png`：三张结果图。
- `coe_reflection_rgb_without_height_phase.png`、`coe_reflection_rgb_with_height_phase.png`：两种高度相位状态的纯 RGB 预览图。
- `coe_reflection_d65_color.png`：大尺寸 D65 合成图，上排全范围、下排局部放大；放大范围由 `colorView` 控制，仅改变显示。
- `coe_reflection_srgb_linear_without_height_phase.png`、`coe_reflection_srgb_linear_with_height_phase.png`：D65 原始线性曝光经 sRGB 编码后的对照图，不使用亮度增强。

远场对比图改为横向排列：上行为“不叠加高度相位”，下行为“叠加高度相位”；每列是一个预览波长，最后一列是光谱色度合成图。单波长图统一灰度 dB 色标和能量参考，色彩图不使用 dB 色标。青线标出雾度半角，黄线标出观察圆锥/环带，两者默认重合；红点线标出零级排除范围，品红圈及 S 编号标出非零级积分域。字体、坐标轴、角度范围和色标样式由脚本统一控制。

D65 模式将逐波长强度与 D65、CIE 三刺激值函数和积分步长相乘后求和，得到 XYZ，再计算 xy 色坐标并转换到 D65 白点的 sRGB。两状态使用相同的曝光和亮度映射，不分别归一化三个通道。图像是当前简化光学模型下的色度预测，不保证等同于真实材料或未校准屏幕的观感。

主脚本提供以下显示参数：

```matlab
colorOptions.fixedGain = 0;        % 共同自动曝光；正数用于固定跨方案曝光
colorOptions.exposureEV = 0;       % 两状态共同曝光补偿
colorOptions.toneMapping = "log";  % log=亮度增强；none=关闭增强
colorOptions.toneStrength = 3000;  % 越大越突出弱衍射
colorOptions.luminanceGamma = 0.75; % 小于 1 提亮暗部，仅作用于亮度 Y
colorView.centerDeg = [0, 0];      % 独立 D65 图下排局部放大的中心
colorView.halfAngleDeg = 1.5;      % 局部显示半范围，必须位于计算窗口内
```

亮度增强按共同曝光后的 Y 计算对数/伽马映射，再以同一个比例缩放该像素的 XYZ。正值高光按整个 RGB 向量等比例压缩，避免逐通道截断引起明显偏色。只要颜色在 sRGB 色域内，这些缩放保持 xy 色度；色域外的负线性 RGB 仍需裁剪，不能保证显示颜色完全不变。显示图是亮度增强预览，不再表达物理亮度的线性比例。

`colorData` 中的 `XYZ`、`xy`、`XYZDensity`、`linearSRGBDensity` 保留原始结果；`sRGB` 为未增强的线性曝光版，`displayXYZ` 和 `displaySRGB` 是单独的显示数据。亮度参数、增强图及缩放不参与任何物理指标计算。

RGB 模式仍用原亮部截断、对数压缩及伽马增强代理预览，参数为 `rgbOptions`。所有显示映射只作用于图片，不改变物理数组。将 `saveResults` 设为 `false` 可只计算和显示。

有限输出窗口可能触发能量不足 98% 的提示。请查看 `CapturedFraction`；如需扩大观察范围，可调整 `angleHalfRangeDeg`（不超过近轴模型的 10° 上限），并结合采样设置检查收敛。累计曲线和观察角域只统计输出窗口内可完整覆盖的角域；全 FFT 总能量仍是近轴模型总量，不代表真实全半球散射。

### 雾度评价值与 5° 圆锥的定义

常用 ASTM D1003 透射雾度的角度阈值是偏离入射方向 **2.5°**，即圆锥半角 2.5°、全角 5°，分母为总透射光。参见 [HunterLab：Measuring Haze](https://support.hunterlab.com/knowledgebase/article/KA-06011/en-us) 和 [ASTM D1003 标准范围](https://store.astm.org/standards/d1003)。

本模块计算的是反射，因此使用 **“雾度评价值（反射锥外能量占比，模型代理指标）”**，而不标为标准透射雾度或标准反射雾度。标准反射雾度使用镜面角附近的特定测量几何，不能直接等同于全角 5° 圆锥之外的全部能量占比，参见 [Rhopoint：Reflection Haze](https://www.rhopointinstruments.com/faqs/appearance-measurement/)。

本模块逐波长的评价公式为 `100*(1 - E_inside_cone/E_total_reflected)`，单位 %。`cfg.metricAngleDeg` 始终表示**半角**，默认 2.5°，圆锥中心是法向入射下的镜面反射方向。分母采用模型全 FFT 网格的反射能量，而不是显示窗口内能量，因此指标包含窗口外能量。它是标量近轴模型的估算，也未做标准仪器所需的光度加权及校正。若没有反射能量，该指标为 NaN（未定义）。

对比表还保存 `TotalReflectancePercent`、`ConeOutsideReflectancePercent`、`ZeroOrderReflectancePercent`，这三项均以**入射能量**为分母，便于区分绝对反射能量和反射能量内部的重新分配。零反射时后两项由未定义的归一化占比计算，仍为 NaN。高度示例会将较多能量移出 ±6° 窗口，`CapturedFraction` 应以实际结果为准。

## 函数接口：供批量计算调用

原有函数接口继续保留。只需添加本方向 `src/` 路径，不添加 `private/`：

```matlab
addpath('topic_coe_reflection/src'); % 在仓库根目录执行；独立复制后改成实际路径
cfg = coeReflectionDefaults();
cfg.featureIntensity = 0.60;
cfg.backgroundIntensity = 0.01;
r = runCOEReflection(cfg);
disp(r.metrics);
imagesc(rad2deg(r.coordinate), rad2deg(r.coordinate), r.kernel(:,:,2));
axis image; colorbar;
```

`heightMap` 以米为单位，反射使用双程相位；支持自定义 `mask` 和复数 `complexCoefficient`。`baselineKernel/baselineMetrics` 为清除额外高度图后的基准。此函数接口只返回物理数组；主脚本负责展示与文件保存。

本目录自带默认参数、FFT/能量积分核心和几何生成代码，不依赖 UDC、旧孔径或根目录辅助函数。法向、标量、近轴远场简化；默认光谱与反射率为示例值。完整指标解释见仓库 `docs/COE反射与UDC成像分支.md`。

## 目录与独立运行

本方向统一使用以下结构：

| 目录 | 用途 |
|---|---|
| `Input/` | 输入图片序列、项目参数文本等计算数据 |
| `Output/` | 计算产生的图片、表格、MAT 等结果 |
| `src/private/` | 辅助函数；保留 MATLAB 特殊目录名，由主程序自动解析 |
| `docs/` | 使用说明、报告等文档 |
| `dist/` | 打包发布文件 |
| `src/` | 主程序源码 |

独立复制整个方向文件夹，在 MATLAB 中仅添加其 `src/` 路径后运行上述入口。输入和输出路径根据入口文件位置定位，与当前工作目录无关。辅助函数目录是命名例外：MATLAB 只有对名为 `private`、位于调用程序同一父目录下的文件夹才会自动应用私有函数解析规则。因此辅助函数放在 `src/private/`，保留小写名称，不改为 `PrivateFunction`，也不单独 `addpath`。

## 按项目切换图片集（2026-11-13）

当前默认使用 `Input/Case1/` 的 PNG 图片集；`Case2` 为 BMP，`Case3` 为 JPG，三者为同几何的演示数据。主脚本第 1 节设置 `projectName`，第 2.2.1 节使用 `stack.geometry.source="images"`。各图层默认不写扩展名，自动匹配该项目文件；同名多格式时必须明确扩展名。BM 白色表示实心材料，`bmImageIsOpening=false`。结果保存至 `Output/项目名称/`。材料 nk 路径也相对于当前项目目录；几何采样和材料厚度仍需用户按项目设置。`projectName=""` 使用旧根目录。详见 `Input/README.md` 和实际输入指南。
