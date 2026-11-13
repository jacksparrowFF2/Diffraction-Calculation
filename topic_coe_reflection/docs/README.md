# COE 反射衍射：使用与交接说明

更新日期：2026-11-13。本文件是本模块的统一使用说明，包含项目配置、实际图样输入、物理模型、结果解释与开发交接。文中 `src/`、`Input/`、`Output/` 均相对于 `topic_coe_reflection/`。

正式仿真读取用户提供的阳极、BM、RGB CF 图样，不调整像素形状。主入口是 `src/COEReflection_main.m`；选定项目后，自动读取该项目的 `config.txt`，无需反复修改 MATLAB 中的常用参数。额外高度默认 **0 nm**。

## 目录

- [快速开始与项目切换](#quick-start)
- [项目配置 config.txt](#project-config)
- [图样、物理尺寸与材料输入](#inputs)
- [光学模型与计算流程](#model)
- [指标与结果解释](#metrics)
- [图像显示与输出文件](#outputs)
- [迁移、验证与开发交接](#handoff)
- [常见问题](#faq)
- [历史示例记录](#history)
- [数据与公式来源](#sources)

<a id="quick-start"></a>

## 快速开始与项目切换

1. 使用 MATLAB 打开 `src/COEReflection_main.m`。本模块已在 Windows / MATLAB R2024a 验证，核心计算不需要 Simulink，也不依赖其它 topic 的源码。
2. 主脚本第 1 节选择 `projectName="Case1"`，保持 `useProjectConfig=true`。
3. 编辑 `Input/Case1/config.txt`，完整运行主脚本。脚本会执行 `clearvars`、`clc`、`close all`，保留搜索路径与 MCP 连接；不要在运行前通过工作区变量设置参数。
4. 结果保存到 `Output/Case1/`，同一 `caseName` 再次运行会覆盖同名文件。

```text
Input/
  Case1/              PNG 示例项目
    anode.png
    BM.png
    CF_R.png
    CF_G.png
    CF_B.png
    config.txt
    README.md         本项目图层含义和采样清单
  Case2/              同几何 BMP 示例
  Case3/              同几何 JPG 示例
  Colorimetry/        共享 CIE 数据，不必复制到每个项目
Output/
  Case1/              Case1 的图像、MAT、CSV 和配置记录
src/
  COEReflection_main.m
  coeReflectionDefaults.m
  runCOEReflection.m
  private/            内部算法和绘图函数
```

三套示例均为 **1024×1024**，Pitch=50 μm，16×16 个 R/B 网格间距，每周期64点，物面800×800 μm。参考图的实线 PDL 作为阳极，虚线是 Mask，不当作 BM。R/B 在网格上交错，G 位于四邻中心并交替倾斜；每个2×2晶胞有2R、4G、2B。BM 开口和 CF 按先前约定分别外扩5 μm和7 μm。参考图未标注尺寸，示例只匹配布局及轮廓特征，不代表产品尺寸或实测材料。

新增项目时复制一个 Case 文件夹并改名，替换五张图片及配置，再修改 `projectName` 即可（支持中文名称）。物理尺度、采样和材料由项目配置指定，不由文件夹名称推断；程序不自动缩放或平铺图片。`projectName=""` 使用旧的 Input/Output 根目录，启用配置时仍需在 Input 根目录提供 config.txt。

默认使用 D65 380–780 nm、5 nm 步长（81波长），FFT补零倍数2，输出401×401、角域±6°。`cfg.heightMap=0`、`useHeightExample=false`；两种高度状态应相同，膜层自身的相位始终保留。

<a id="project-config"></a>

## 项目配置 config.txt

配置在光谱、几何和 FFT 计算前加载，**覆盖主脚本默认值**。以下示例应写入项目 config.txt；设置 `useProjectConfig=false` 时才完全采用主脚本参数。完整常用配置已经提供在 Case1~3 中。

常改参数如下（仅示例，完整配置已放入三个 Case）：

```text
# 物理尺度：单位 m。所有图片覆盖完整物面，不自动缩放和平铺。
cfg.pitch = 50e-6
cfg.periods = [16, 16]
cfg.samplesPerPeriod = 64

# 图片文件名可不写扩展名
stack.geometry.anodeFile = "anode"
stack.geometry.bmFile = "BM"
stack.geometry.bmImageIsOpening = false
stack.geometry.cfFiles = ["CF_R", "CF_G", "CF_B"]

# 材料 nk 和厚度；文件非空时优先于 n/k 常数
stack.materials.BM.nkFile = "Materials/BM_nk.csv"
stack.materials.BM.thickness = 1e-6

# 光谱、观察角域与显示
spectrumMode = "D65"
wavelengthStepNm = 5
observation.halfAngleDeg = 2.5
colorOptions.exposureEV = 0
colorOptions.toneStrength = 3000
displayStyle.showAnnotations = true
cfg.heightMap = 0
```

### 文本规则

- UTF-8 编码，每行 `字段名 = 值`。空行及以 `#`、`%` 开头的整行注释会跳过。
- 数值使用十进制或科学计数，如 `50e-6`；长度单位 m，角度单位 °，光谱步长及预览波长单位 nm。
- 数组采用 JSON 格式，英文逗号分隔，例如 `[16,16]`；布尔值为 `true`、`false`；字符串和路径必须用双引号。Windows 子路径推荐使用 `/`。
- 不支持行尾注释、分号、`50*1e-6`、`Inf`、`eye(3)` 或其它 MATLAB 表达式；不会执行文本里的代码。
- 可只保留需要覆盖的字段；其余采用主脚本默认值。只支持已有参数，不允许添加未知字段，也不允许配置像素尺寸、圆角和倾角。
- 改预览波长时同时检查 `spotOptions.wavelengthNm`；该参数独立控制光斑统计波长。
- 采样参数需与输入图片匹配。切换项目不改变图片本身；程序不会根据配置重画图样。

结果 MAT 保存最终 cfg、stack 以及加载的配置文本，输出目录还会保存 `*_config.txt`，便于回溯。同名项目输出仍会覆盖。共享 CIE 数据路径不受项目切换影响。

需临时完全采用主脚本参数时设置 `useProjectConfig=false`；此时不读取 config.txt。默认应保持 true。

<a id="inputs"></a>

## 图样、物理尺寸与材料输入

所有图样和 nk 路径均相对于当前 `Input/项目名称/`，可按需使用 `Layers/`、`Materials/` 子目录。以下列出的参数均可写入 config.txt；主脚本内同名设置只作为默认值。

### 先确定图样物理尺寸与采样

所有层必须覆盖同一个完整矩形物面、采用相同原点与像素网格，行为 y、列为 x。不自动识别 Pitch，不缩放图片，也不把单个像素或晶胞自动平铺。需要重复版图时，先在外部生成完整矩阵。

以下是对应的 MATLAB 参数写法（复制到 config.txt 时去掉分号和行尾注释）：

```matlab
cfg.pitch = 50e-6;            % 相邻像素单元间距 50 μm；不是 RGB 晶胞的总尺寸
cfg.periods = [16, 16];       % [x 方向单元数, y 方向单元数]
cfg.samplesPerPeriod = 64;    % 每单元的等间距采样点数
```

对应所有图样均为 **1024 行 × 1024 列**，采样间距 **0.78125 μm**，物面总宽/高为 **800 μm**。通式：

- 列数 Nx = `cfg.periods(1)*cfg.samplesPerPeriod`。
- 行数 Ny = `cfg.periods(2)*cfg.samplesPerPeriod`。
- 采样间距 dx = dy = `cfg.pitch/cfg.samplesPerPeriod`。
- 物面尺寸 Lx = `cfg.periods(1)*cfg.pitch`，Ly = `cfg.periods(2)*cfg.pitch`。

当前使用正方形采样像素。`periods` 采用偶数与现有配置保持一致。若导出的图样每个 50 μm 单元有 100 点，则将 `samplesPerPeriod` 改为 100，并按上述公式核对完整矩阵尺寸。需要保留小于原采样间距的窄边/间隙时，应先提高真实版图的栅格采样；`paddingFactor` 只加密远场采样，不能修复物面欠采样。

### 推荐：一个 MAT 文件导入各层

配置保持 `inputFile=""`、`cfg.mask=[]`、`cfg.complexCoefficient=[]`，并设置：

```matlab
opticalModel = "multilayer";
stack.geometry.source = "mat";
stack.geometry.matFile = "Layers/coe_layers.mat";
cfg.heightMap = 0;
useHeightExample = false;
```

MAT 文件顶层字段如下，所有矩阵均为 Ny×Nx 的 `logical` 或严格数值 0/1：

| 字段 | 是否必需 | 数值 1 的含义 |
|---|---|---|
| `anodeMask` | 是 | 存在阳极反射底 |
| `cfRMask` | 是 | 存在 R CF 膜 |
| `cfGMask` | 是 | 存在 G CF 膜 |
| `cfBMask` | 是 | 存在 B CF 膜 |
| `bmOpenMask` | 可选 | BM 开口，没有 BM 材料 |
| `bmMask` | 可选 | 实心 BM 材料 |

推荐在导出程序中保存：

```matlab
% 各变量已由真实版图生成；不在本段凭空创建几何。
save('coe_layers.mat', 'anodeMask', 'bmOpenMask', ...
    'cfRMask', 'cfGMask', 'cfBMask');
```

`bmOpenMask`、`bmMask` 任选其一；同时提供时必须互为补集。均未提供时，才按 `stack.bmExpansion` 对阳极外扩生成 BM 开口。因此**导入真实 BM 图样后，修改 bmExpansion 不会再改变该图样**。`bmImageIsOpening` 只用于图片，不影响 MAT 字段含义。

### 或者：五张黑白图片

支持 PNG、BMP、JPG/JPEG、TIF/TIFF。推荐无透明度、无调色板的 8 位灰度 PNG 或 BMP；二值 TIFF 也可。勿用带坐标轴、标尺、文字的截图。JPG 可以读取，但有损压缩可能影响几何边缘，实际二值图样优先用 PNG/BMP。各图大小必须正好为 Ny×Nx。程序用 `imread` 读取、归一化后，以严格大于 `threshold` 的位置为 1；RGB 图会先对三个通道取平均，所以**不要直接输入一张红绿蓝合并图**，分别输出每个 CF 的黑白掩膜。

```matlab
stack.geometry.source = "images";
stack.geometry.anodeFile = "anode";  % 白色=阳极
stack.geometry.bmFile = "BM";
stack.geometry.bmImageIsOpening = false;      % 白色=实心 BM；白色=开口时改 true
stack.geometry.cfFiles = ["CF_R", "CF_G", "CF_B"];
stack.geometry.threshold = 0.5;
```

不带扩展名时，程序自动查找项目内的对应图片；同名多格式会要求明确扩展名。也可直接指定如 `anode.bmp`。每套示例都有五张 1024×1024 图片，只需修改项目名称即可体验格式切换。

CF 图片白色表示有该色膜，黑色表示无该色膜。MAT 和图片输入均不翻转坐标；核对输出的各层几何图与坐标方向。若不提供 BM 图片，将 `bmFile=""`，使用 `bmExpansion` 自动生成。

当前不直接导入 DXF、DWG、GDS 等矢量版图，应先按统一物理尺度和原点栅格化为上述 MAT 或图片。

### 材料 nk 数据与厚度

每种材料一个 CSV，推荐**无表头、无单位行、无空行**，恰好三列：`波长 nm,n,k`。例如文件的格式为：

```csv
380,1.72,0.08
400,1.70,0.06
550,1.65,0.03
780,1.61,0.02
```

这些数值仅说明格式，请替换为实测数据。波长严格递增、不重复，n>0、k≥0、数值有限。使用 `N=n+i*k` 的约定；不要把吸收系数或介电常数直接放进 k 列。波长必须以 **nm** 存储；数据须覆盖全部计算波段（默认 380–780 nm），不作外推。程序将 n 和 k 分别线性插值至 `cfg.wavelengths`，不要求原始数据等间距，也不要求与 5 nm 计算网格一致。

对应的 MATLAB 参数写法如下（config.txt 中需去掉分号和行尾注释）：

```matlab
stack.materials.BM.nkFile = "Materials/BM_nk.csv";
stack.materials.R.nkFile  = "Materials/CF_R_nk.csv";
stack.materials.G.nkFile  = "Materials/CF_G_nk.csv";
stack.materials.B.nkFile  = "Materials/CF_B_nk.csv";

% 以下厚度数值只是填写示例，请换成各层实际厚度；单位 m。
stack.materials.BM.thickness = 1000e-9;
stack.materials.R.thickness  = 1100e-9;
stack.materials.G.thickness  = 1000e-9;
stack.materials.B.thickness  = 900e-9;
```

`nkFile` 非空时优先于该材料的 `n`、`k` 常数/向量；厚度仍由 `thickness` 单独指定，不写在 CSV 第四列。`nkFile=""` 则使用脚本中的 n/k，可为常数或与计算波长一一对应的向量。

阳极当前不输入 nk 或厚度：保留 `stack.anodeReflectance=1`、`stack.anodePhaseDeg=180`。阳极外是匹配吸收底，绝不扩大阳极反射范围；上覆 BM/CF 的界面仍会产生反射。

### 检查和输出

自动检查 BM 开口包含阳极、RGB CF 不交叠，且默认 `requireCFCoverage=true` 时 CF 完整覆盖 BM 开口、每个非空 CF 向 BM 实体上方有所延伸。不要为绕过真实版图对齐错误而关闭检查。

运行前可设 `caseName="actual_coe_01"` 区分输出；同名文件会覆盖。优先核对输出几何图、`*_layer_patterns.mat`、`*_stack_materials.csv`，确认最终掩膜、插值后的 n/k 和厚度。然后查看 D65 图、累计角能量与绝对反射率 CSV。显示增强不改变原始 XYZ 和物理指标。

额外高度保持 0 nm 时，两种高度状态应一致；膜层自身的厚度干涉相位始终保留。当前仍是法向局部相干薄膜 + 标量近轴远场模型，忽略侧壁和横向电磁耦合。

<a id="model"></a>

## 光学模型与计算流程

主脚本按“初始化 → 读取项目配置 → 读取图样/材料 → 计算叠加高度状态 → 清零额外高度并再次计算 → 积分评价 → 展示/保存”顺序组织。具体算法在 private 中，用户入口保持完整、可从上往下阅读。

### 阳极、BM 与 CF

- 平面像素采用 Diamond 排布；用户提供 Pitch，例如 50 μm。
- 从下到上为阳极、BM、RGB CF。各层是图案化区域，并非整面均匀覆盖。
- **阳极反射只存在于输入阳极图样范围内**，其裸底反射率默认 100%。
- BM 开口由阳极开口边界向外扩张 5 μm 得到，外扩量可调；开口之外是实心 BM。
- **BM 开口中超出阳极图样的外扩环带属于不透明吸光底层，不能按 100% 阳极处理。**
- R/G/B CF 覆盖各自对应的 BM 开口，其区域大于 BM 开口；各色 CF 区域之间不得交叠。
- 用户可输入阳极、BM 和 RGB CF 图样，并提供 BM / R CF / G CF / B CF 的 n、k 与厚度。

### 局部膜堆

| 区域 | 从入射侧向下的局部结构 |
|---|---|
| 阳极内、BM 开口内、有 CF | 空气 / 对应 CF / 阳极反射底 |
| 阳极外、BM 外扩开口内、有 CF | 空气 / 对应 CF / 吸光底层 |
| BM 实体上方有 CF | 空气 / 对应 CF / BM / 该位置的底层 |
| BM 实体上方没有 CF | 空气 / BM / 该位置的底层 |

底层由阳极掩膜判定，不能仅依据 BM 开口判定。阳极默认裸底 r=-1（R=100%、相位180°），反射率和相位均可调。阳极外的理想吸收终端在直接接触入射介质时裸底 r=0，即阻抗匹配；进入终端的净能量被吸收，不作为背后透射光再次返回。加上 CF/BM 后，按同一个终端阻抗重新计算膜层/终端界面的反射，所以局部最终系数不一定为零。这是理想边界，不是已经测得的吸光材料 n/k。

输入图样必须有一致的物理尺度、坐标原点和矩阵尺寸。分别保留阳极掩膜、BM 开口掩膜和三个 CF 掩膜；检查 BM 包含阳极、CF 覆盖开口及 CF 互不交叠。若导入实际 BM 图样，不应静默用自动外扩结果替换它。

各位置先按局部平行平面膜堆计算随波长变化的复反射系数，再进入已有 FFT / 角域 / CIE 后处理。CF 局部厚度固定，位于 BM 实体上方时叠加在 BM 顶面，在开口内则从裸底平面起算。以 `BM厚度+最大CF厚度` 为共同空气参考面，补足各柱体顶面到该平面的往返相位。未包含侧壁或邻近柱体的矢量耦合，不能当作 RCWA/FDTD。

膜厚传播与干涉相位属于膜层模型，不等同于可选的额外高度图；`heightMap` 默认仍为 0。几何采用 regionMap，局部光谱系数采用 16×波长数的查找表，避免构建庞大的逐像素光谱系数三维数组。



该模型为法向局部平行平面相干膜堆 + 标量近轴远场，忽略侧壁、邻近柱体矢量耦合、斜入射、偏振、有限角尺寸光源等因素，不是 RCWA/FDTD 或真实全半球散射。常数 CF n/k 不保证实际 RGB 通带，需输入真实色散吸收数据。

### 光谱与色度积分

令 S_i 为 D65 相对光谱、q_i 为梯形积分宽度（nm）。各波长的输出数组已相对该波长入射能量归一化。

- 辐射度：`w_i = S_i*q_i / sum(S*q)`。
- 明视觉：`wY_i = S_i*q_i*ybar_i / sum(S*q*ybar)`。
- XYZ：`WXYZ_i = S_i*q_i*[xbar_i,ybar_i,zbar_i] / sum(S*q*ybar)`，三个分量共用一个归一化常数。

因此辐射度和明视觉绝对反射率分别以入射辐射能与入射明视觉光通量为分母。二者不能混为一项，也不能先对每个波长的“反射能量内部比例”求平均来替代光谱能量积分。

所有光谱量只覆盖 380–780 nm，未把波段外的入射能量计入分母。D65 在该范围、5 nm 步长下得到白点 xy 约 (0.31272, 0.32904)，与 D65 标准色度接近；这是数据及积分的检查，不代表波长步长已经对所有尖锐衍射结构收敛。

### 两种高度状态与兼容入口

“不叠加高度相位”仅将额外 heightMap 清零，“叠加高度相位”采用输入的 heightMap。两者仍保留膜厚干涉和共同参考面的相位；它们是同一输入结构的两种计算状态。`baseline` 字段名为兼容历史接口保留。

额外高度相位 `phi=4*pi*deltaIndex*h/lambda`，h单位m。均匀高度只增加整体相位，非均匀高度才重新分配远场能量；`deltaIndex` 是额外光程层的有效折射率差，不能直接当作某一种膜的n。不要把膜厚相位重复加入 heightMap。

可选演示 `useHeightExample=true` 使用0/687.5 nm二值条纹，默认周期Pitch/4=12.5 μm、占空比0.5、沿x变化，550 nm下相位差为π（Δn=0.20）。它不是当前默认值或产品优化方案；周期至少4个采样点，周期和高区宽度须与采样间距整除。非空 inputFile 时跳过演示。

`inputFile` 是旧单层/额外高度/外部复系数 MAT 入口，与多层图样的 `stack.geometry.matFile` 不同。可包含 mask、heightMap、complexCoefficient：mask为二维0–1区域图，高度为标量或完整二维矩阵，复系数可为标量、二维或Ny×Nx×波长数，模不大于1。多层模式请通过 stack.geometry 输入各层，不使用旧 mask。外部复系数非空时优先于多层模型；若其已包含额外高度相位，应设置 heightMap=0。

`opticalModel="simple"` 保留旧单界面模型。启用 cfg.materials 时使用 `r=(n_incident-(n+i*k))/(n_incident+(n+i*k))`；关闭时使用 featureIntensity/backgroundIntensity 的强度反射率（振幅取平方根）。`spectrumMode="RGB"` 保留三波长代理预览，不是标准色度/光度计算。正式项目仍采用用户图样，不提供像素形状配置。

<a id="metrics"></a>

## 指标与结果解释

`spectralEnergy(:,:,k)` 是每个输出角度 bin 的积分能量，相对于该波长入射能量；它不是未积分的强度密度，也不能直接命名为 BRDF。物理指标均来自原始数组，不使用显示增强图。

| 指标 | 分母与用途 |
|---|---|
| 总反射率 Throughput | 模型总反射能量/入射能量，判断反射总量 |
| 观察角域绝对反射率 | 视锥或环带内反射能量/入射能量；明视觉版本按相应光通量定义 |
| 累计角能量 | 指定中心的实心圆锥积分/入射量，描述能量如何随半角累计 |
| 非零级光斑绝对反射率 | 每个固定积分圆内能量/入射量，同时报告角位置 |
| 零级与锥外占比 | 各自能量/模型总反射能量；要与绝对反射率一起解释 |
| CapturedFraction | 显示窗口内能量/模型总反射能量，检查有限窗口损失 |
| ParsevalAbsoluteError | FFT 能量守恒绝对误差 |

比例字段一般为0–1，名称带 Percent 的汇总列为百分数。零反射能量时归一化占比为 NaN，不代表零雾度。纯相位可能降低零级但增强其它方向亮斑，因此不能只优化零级保持率或锥外占比。

### 观察角域与累计曲线

`observation.centerDeg` 定义 [theta_x, theta_y] 中心，`halfAngleDeg` 为外半角，`innerHalfAngleDeg` 为内半角；内半角为 0 表示实心圆锥，正值表示环带。圆锥中心可偏离镜面方向。要求整个角域位于输出窗口内，越界会报错而不是静默丢失一部分能量。

每个波长先在指定角域积分，再按上节权重合成宽带结果。累计曲线则是围绕该中心、由 0 增至输出窗口内最大完整圆半径的实心圆锥积分，分母始终为入射量。零半径定义为 0；其它半径按输出 bin 中心作边界判断，应检查空间采样及角度采样收敛。

对于观察环带，表格直接统计环带能量，而累计曲线仍是从中心起算的实心圆积分。曲线末端既不包含方形输出窗口的所有角点，也不包含窗口之外的能量，不能将其强制归一化为 100%。

### 非零级衍射斑

主脚本的 `spotOptions` 可设置零级排除半角（默认 0.30°）、积分半角（0.15°）、最小中心间距（0.35°）、每组最多光斑数（8）和相对峰阈值（0.08）。`wavelengthNm` 默认统计预览波长，空数组可统计全部波长，另始终计算宽带 `radiant/photopic` 两组。

检测基于两状态原始能量图的共同候选，分别在同一组圆域中积分，输出以入射量为分母的绝对反射率。圆域不进入零级排除范围、不重叠、不越过输出窗口。图中品红圈和 S 编号对应 `_diffraction_spots.csv`。编号不是理论衍射阶次，不同组的能量不能直接相加。

积分半角不得小于输出角网格间距，最小中心间距至少两倍积分半角。候选积分圆不得跨入零级排除区、相互重叠或越出窗口。宽带色散带不能由一个局部圆代表全部能量；不同组的 S 编号不表示同一理论阶，对称峰编号也可能随浮点误差互换。

### 雾度代理的定义与限制

常用 ASTM D1003 透射雾度的角度阈值是偏离入射方向 **2.5°**，即圆锥半角 2.5°、全角 5°，分母为总透射光。参见 [HunterLab：Measuring Haze](https://support.hunterlab.com/knowledgebase/article/KA-06011/en-us) 和 [ASTM D1003 标准范围](https://store.astm.org/standards/d1003)。

本模块计算的是反射，因此使用 **“雾度评价值（反射锥外能量占比，模型代理指标）”**，而不标为标准透射雾度或标准反射雾度。标准反射雾度使用镜面角附近的特定测量几何，不能直接等同于全角 5° 圆锥之外的全部能量占比，参见 [Rhopoint：Reflection Haze](https://www.rhopointinstruments.com/faqs/appearance-measurement/)。

本模块逐波长的评价公式为 `100*(1 - E_inside_cone/E_total_reflected)`，单位 %。`cfg.metricAngleDeg` 始终表示**半角**，默认 2.5°，圆锥中心是法向入射下的镜面反射方向。分母采用模型全 FFT 网格的反射能量，而不是显示窗口内能量，因此指标包含窗口外能量。它是标量近轴模型的估算，也未做标准仪器所需的光度加权及校正。若没有反射能量，该指标为 NaN（未定义）。

对比表还保存 `TotalReflectancePercent`、`ConeOutsideReflectancePercent`、`ZeroOrderReflectancePercent`，这三项均以**入射能量**为分母，便于区分绝对反射能量和反射能量内部的重新分配。零反射时后两项由未定义的归一化占比计算，仍为 NaN。高度示例会将较多能量移出 ±6° 窗口，`CapturedFraction` 应以实际结果为准。

<a id="outputs"></a>

## 图像显示与输出文件

远场对比图横向排列，上行为“不叠加高度相位”，下行为“叠加高度相位”；列为预览波长与 D65 合成图。D65 专用图包含全范围和局部放大，后者仅裁切显示，不修改计算和积分范围。单波长统一灰度dB色标与参考值，色彩图不使用dB色标。

### 共同曝光与亮度增强

对每个输出角度 bin，按 XYZ 权重合成三刺激值，计算 `x=X/(X+Y+Z)`、`y=Y/(X+Y+Z)`。无信号处 xy 为 NaN。保存的 XYZ 是该 bin 的积分量，以单位入射明视觉量为参考。

显示前将 XYZ 除以近轴角度 bin 面积，得到每单位近轴立体角的信号密度，避免输出网格变密导致图像自动变暗。这是归一化角分布，不直接等同于实测 BRDF 或 cd/m²。

XYZ 采用 D65 白点的 XYZ→线性 sRGB 矩阵转换，再使用同一个曝光倍数和标准 sRGB 分段编码。`colorOptions.fixedGain=0` 时，以两状态共同的最大线性 RGB 密度为参考；设为正数时固定曝光，适合不同设计间保持一致。`exposureEV` 对两状态共同乘以 2^EV。

现在默认启用可调的亮度增强，不做逐通道增强、亮度匹配或白平衡。对共同曝光后的亮度 Y 使用 `Y_target=[log(1+alpha*Y)/log(1+alpha)]^gamma`，再对 XYZ 统一乘以 `Y_target/Y`。默认 alpha=3000、gamma=0.75；两种状态使用同一个函数和同一组参数，设 `toneMapping="none"` 可关闭。

高光通过 RGB 三个分量统一除以 `max(1,max(R,G,B))` 压入可显示范围，XYZ 同步缩放。在 sRGB 色域内，以上处理保持 xy 不变；色域外负线性 RGB 的裁剪仍可能改变显示色度。增强显示不保留亮度的线性比例，不应拿增强图像的亮度计算能量指标。

原始 `XYZ`、`xy`、`XYZDensity`、`linearSRGBDensity` 均不修改，`sRGB` 保留原线性曝光图；`displayXYZ` 和 `displaySRGB` 单独保存增强结果。`gamutClippedFraction` 记录原线性曝光的裁剪像素比例，`displayGamutClippedFraction` 记录增强及整体高光压缩后仍需裁剪的像素比例，二者都不是能量损失比例。若要精确比较颜色，应基于保存的原始 XYZ/xy，而不是增强截图。

`coe_reflection_d65_color.png` 另外提供两状态的全范围和局部放大图。默认局部显示 ±1.5°，可通过 `colorView.centerDeg` 和 `colorView.halfAngleDeg` 调整；不改变仿真采样、物理观察角域或累计角能量统计。

该图是基于 CIE 的色度预测，不再是三条窄带线直接充当 R/G/B 的代理颜色；但材料常数、简化光学模型、未知绝对照明和显示设备仍限制它与真实产品观感的一致性。

显示参数包括 `colorOptions.fixedGain`（0为共同自动曝光，正数固定曝光）、`exposureEV`、`toneMapping`（log或none）、`toneStrength`（默认3000）、`luminanceGamma`（默认0.75），以及 `colorView.centerDeg/halfAngleDeg`（默认中心0、局部±1.5°）。跨项目比较应固定或记录曝光。RGB代理模式另使用 rgbOptions 的亮部截断、对数和伽马参数。

### 图形标记开关

标记在绘图阶段叠加。`displayStyle.showAnnotations=false` 隐藏全部辅助标记；细分开关为 showMetricCone（青色评价锥）、showObservationRegion（黄色观察角域）、showZeroOrderExclusion（红色零级排除）、showSpotRegions（品红积分圆）、showSpotLabels（S编号与引线）。这些开关只影响图形，不改变检测、积分、原始XYZ和CSV。

### 保存文件

以下后缀均加上 caseName 前缀，保存在 Output/项目名称/。`saveResults=false` 只计算和展示。
| 文件后缀 | 内容 |
|---|---|
| `_results.mat` | 实际参数、原始光谱能量、XYZ/xy、增强图、所有统计结构 |
| `_diffraction.png` | 横向逐波长 / D65 对比图，含光斑标记 |
| `_d65_color.png` | 大尺寸 D65 全范围与局部图，含光斑标记 |
| `_geometry.png`、`_metrics.png` | 几何及零级/雾度代理指标 |
| `_cumulative_energy.png/.csv` | 累计绝对反射率 |
| `_observation_metrics.csv` | 观察角域绝对反射率和该区域色坐标 |
| `_diffraction_spots.csv` | 非零级光斑中心、实际峰位置、积分半角、绝对反射率 |
| `_metrics.csv`、`_baseline_metrics.csv` | 两状态的逐波长物理指标 |
| `_comparison_metrics.csv` | 两状态的角域、零级及锥外指标对比 |
| `_spectrum.csv` | D65 与三刺激值函数的实际采样数据 |
| `_stack_regions.csv`、`_stack_materials.csv` | 多层区域编码、阳极/BM/CF 状态与面积比例；实际逐波长材料 n/k、厚度 |
| `_layer_patterns.mat` | 实际阳极、BM 开口和 RGB CF 图样，可复制到 Input 复用 |
| `_rgb_without/with_height_phase.png` | 纯增强 RGB 图，无标记 |
| `_srgb_linear_without/with_height_phase.png` | 未增强的线性曝光版 sRGB 图，无标记 |

光斑表中 `Group="wavelength"` 的分母为该波长入射能量；`radiant` 为 D65 波段总入射辐射能；`photopic` 为 D65 入射明视觉光通量（RGB 模式是代理权重）。`AbsoluteReflectancePercent=100*E_ROI/E_incident`，不是光斑在总反射能量中的占比。宽带行 `WavelengthNm=NaN` 是正常的，表示多波长合成。不同组、不同波长、两个状态之间**不能相加**。

`*_config.txt` 保存本次加载的配置文本；`*_results.mat` 同时保存最终 cfg、stack、projectConfig 以及原始数组。无标记的 RGB/sRGB 图片与增强图均不能用于能量计算。

<a id="handoff"></a>

## 迁移、验证与开发交接

### 在另一台电脑继续工作

复制完整 topic_coe_reflection 文件夹或从 Git 拉取 main，保留 src/private、项目输入和 Input/Colorimetry 的层级即可，不需要迁移 ZIP、MATLAB 安装或原电脑账户配置。输出目录默认忽略于 Git；需要相同原始结果时另行复制 MAT，或在新电脑重跑。

仅将 src 加入 MATLAB 路径，不要单独 addpath(private) 或用 genpath 展开私有目录。路径按主脚本自身位置解析，与工作目录和盘符无关。Codex 控制现有 MATLAB 时需配置新电脑的 MCP，并执行 shareMATLABSession()；手动运行程序不需要 MCP。

### 验证

完整仓库的 COEMultilayerTest、COESpectralTest、COESpotsTest 覆盖膜堆解析极限、独立特征矩阵、吸收/守恒、图样与nk读取、色度、角域和光斑；ModuleIsolationTest 的 COE 分支验证独立复制后的主脚本与输出。

在仓库根目录执行：

```matlab
r = runtests({'tests/COEMultilayerTest.m', ...
    'tests/COESpectralTest.m', 'tests/COESpotsTest.m'});
assert(all([r.Passed]));
suite = matlab.unittest.TestSuite.fromFile('tests/ModuleIsolationTest.m');
suite = suite(contains({suite.Name}, '/coe', 'IgnoreCase', true));
r = run(suite);
assert(all([r.Passed]));
```

默认0 nm高度时两种状态一致；D65白点约(0.31272,0.32904)。验收应核对输入掩膜、最终材料表、能量守恒和当前配置，不把历史几何的反射率当作新版示例固定基准。对称峰按角位置和积分域匹配，不以S编号作为永久标识。最终产品还需检查物面采样、补零、输出角网格和光谱步长收敛并与实验对比。

### 程序责任与后续维护

| 文件或目录 | 用途 |
|---|---|
| src/COEReflection_main.m | 输入、配置、两状态计算、展示和保存 |
| src/coeReflectionDefaults.m、src/runCOEReflection.m | 保留兼容的公共配置/批量计算接口 |
| src/private/loadCoeProjectConfig.m | 项目文本配置解析 |
| src/private/prepareCoeStack.m、reflectCoeFilmStack.m | 图样校验、区域查找表和相干膜堆 |
| src/private/simulateDisplayDiffraction.m | 复振幅、FFT、保能量重采样与基础指标 |
| src/private/loadCoeSpectrum.m、synthesizeCoeColor.m | D65、CIE积分和色度显示 |
| src/private/computeCoeAngularMetrics.m、computeCoeDiffractionSpots.m | 角域与光斑积分 |
| src/private/format*.m、markDiffractionSpots.m、styleCoeFigure.m | 图形样式和辅助标记 |
| Input/Colorimetry/ | 原始CIE数据、元数据、许可与校验值 |
| AGENTS.md | 本模块维护规则 |
| dist/ | 历史打包工具与旧快照，继续工作以当前源码为准 |

src 顶层仅保留主脚本与两个兼容公共接口。新增内部函数直接放在 private，使用小驼峰命名，文件名与首函数名一致；用户主流程和必要中文注释留在主脚本。现有数值算法回归由 tests/private/addCoeAlgorithmFixture.m 在临时副本运行，生产目录的私有解析由独立模块测试覆盖。

继续开发时优先接入实测图样与色散数据，核对CF剖面/平坦化和吸光底参数；可再增加手动ROI/阶次跟踪、角度截线与色度分布、有限光源/入射角/偏振/接收孔径。标准反射雾度需单独实现相应测量几何与校正。用户说明统一维护在本 README，不再按功能拆出多个平行说明文件。

批量接口示例（仅返回物理数组，不自动执行主脚本全部后处理，也不自动加载项目config）：

```matlab
addpath('topic_coe_reflection/src');
cfg = coeReflectionDefaults();
cfg.complexCoefficient = 0.5; % 已知输入复振幅的演示，非实际材料参数
r = runCOEReflection(cfg);
disp(r.metrics);
```

<a id="faq"></a>

## 常见问题

- 默认两个状态相同：这是 0 nm 高度的预期行为，不是计算失效。
- 零级降低但光晕增强：纯相位会重新分配能量，不一定降低总反射率；应结合观察角域、累计曲线和非零级能量判断。
- 找不到光斑：检查阈值、零级排除半径和积分半角，确认候选完整位于输出范围；没有候选时会输出只有表头的 CSV。
- CSV 中没有某些弱光斑：受到相对峰阈值、最小间距、最多光斑数及有限窗口筛选。这些选定光斑不代表所有非零级反射能量。
- 窗口能量不足 98%：查看 `CapturedFraction`，扩大角度范围并检查采样；近轴模型限定半范围不超过 10°。显示窗口外能量不会出现在光斑统计中。
- MATLAB 旧版本函数不支持：以已验证的 R2024a 环境为基准，不直接删除输入检查；先核对版本及语法兼容性。
- 中文字体异常：修改 `displayStyle.fontName` 为新电脑已安装的中文字体。
- 修改后结果异常：先关闭高度、保留默认参数运行验证；检查 m/nm、全角/半角、矩阵行列和波长顺序。

- 改主脚本参数却不生效：当前config会覆盖默认值，请改项目配置，或显式关闭useProjectConfig。
- 图片找不到或同名多格式：核对项目目录和文件名，重复格式时明确扩展名；所有层须同尺寸同原点。
- BM/CF校验失败：核对白色含义、BM开口包含阳极、CF覆盖开口且各色不交叠，不要直接关闭检查掩盖对齐错误。
- 字段/数值报错：核对config的JSON字面量语法、大小写、数字单位及nk波段覆盖范围。

<a id="history"></a>

## 历史示例记录

以下仅保留早期运行记录，**不作为当前 Diamond 图片项目的验收值或完整复现配置**。旧文档曾随默认采样更新而混入不同版本的尺寸描述，因此此处不将其声明为当前参数。精确复现须从 Git 历史取出对应源码、几何和全部配置，不能仅更换光谱选项。

### 三波长高度条纹演示（2026-11-12）

旧simple模型使用0.60/0.01强度反射率、Δn=0.20，512×512物面、补零4倍、401×401输出、±6°窗口；高度为0/687.5 nm、周期12.5 μm、占空比0.5。表中箭头为不叠加→叠加额外高度，零级和锥外均以总反射能量为分母。
| 波长 | 零级能量占比 | 2.5° 半角圆锥外能量占比 | ±6° 窗口捕获比例 |
|---|---|---|---|
| 620 nm | 39.734 → 1.2534 | 11.017 → 84.487 | 97.033% → 83.973% |
| 550 nm | 39.756 → 0.017302 | 9.8151 → 79.837 | 97.768% → 84.254% |
| 450 nm | 39.734 → 4.6632 | 6.6734 → 33.562 | 98.589% → 89.345% |

旧示例两状态总反射率约8.4326%，说明纯相位重分配能量而不减少总反射。一级偏转约λ/周期，RGB分别为±2.84°、±2.52°、±2.06°，550 nm一级靠近2.5°阈值，对角域边界与采样敏感。

### 旧D65单界面高度演示

| 指标，单位 % | 不叠加高度相位 | 叠加高度相位 |
|---|---:|---:|
| 观察角域绝对反射率，辐射度 | 4.50313 | 2.03050 |
| 观察角域绝对反射率，明视觉 | 4.50282 | 1.64783 |
| 模型总反射率，辐射度 | 4.52816 | 4.52816 |
| 模型总反射率，明视觉 | 4.52816 | 4.52816 |

旧多层示例曾记录总辐射/明视觉反射率约9.93634%/9.89988%，中心2.5°内约8.78640%/8.63085%；这些属于修正Diamond布局之前的几何，不能用作当前图片项目基准。

<a id="sources"></a>

## 数据与公式来源

CIE原始CSV、元数据、许可和MD5清单保留在 [Input/Colorimetry/README.md](../Input/Colorimetry/README.md)，不要改写原始数据。CIE数据作者/出版者为International Commission on Illumination，2019，数据遵循CC BY-SA 4.0；此许可不代表其它代码的许可。

- [CIE D65 数据](https://doi.org/10.25039/CIE.DS.hjfjmt59)。
- [CIE 1931 2°标准观察者](https://doi.org/10.25039/CIE.DS.xvudnb9b)。
- [ICC sRGB 定义](https://registry.color.org/rgb-registry/srgb)。
- [Byrnes：Multilayer optical calculations](https://arxiv.org/abs/1603.02720)。
- [MIT：Fresnel interface coefficients](https://web.mit.edu/2.710/OldFiles/Fall06/2.710-wk7-b-sl.pdf)。
- [ASTM D1003 范围](https://store.astm.org/standards/d1003)、[HunterLab：Measuring Haze](https://support.hunterlab.com/knowledgebase/article/KA-06011/en-us)。
- [Rhopoint：Reflection Haze](https://www.rhopointinstruments.com/faqs/appearance-measurement/)、[Rhopoint IQ 使用说明](https://www.rhopointinstruments.com/wp-content/uploads/2023/10/rhopoint-iq-manual-english.pdf)。
- [NIST：Display Reflectance Model Based on the BRDF](https://www.nist.gov/publications/display-reflectance-model-based-brdf)、[Computation of Light Scattered Into Detector](https://www.nist.gov/publications/computation-light-scattered-detector)。
