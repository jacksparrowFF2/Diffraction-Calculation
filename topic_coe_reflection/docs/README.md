# COE 反射衍射：使用与交接说明

更新日期：2026-11-13。本文件是本模块的统一使用说明，包含项目配置、实际图样输入、物理模型、结果解释与开发交接。文中 `src/`、`Input/`、`Output/` 均相对于 `topic_coe_reflection/`。

正式仿真按图层列表读取用户提供的任意命名图样（默认阳极、BM、RGB CF），不调整像素形状。主入口是 `src/COEReflection_main.m`；选定项目后，自动读取该项目的 `config.txt`，无需反复修改 MATLAB 中的常用参数。额外高度默认 **0 nm**。

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
    IJP.png          全白、厚10 μm、n=1.5、k=0
    BM.png
    CF_R.png
    CF_G.png
    CF_B.png
    config.txt
    README.md         本项目图层含义和采样清单
  Case2/              同几何 BMP 示例
  Case3/              同几何 JPG 示例
  CaseSingleNK/       空气/有限膜/空气，只有 singlelayer.bmp
  CaseSingleReflectance/  直接反射率，同样只读一张图
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

默认使用 D65 380–780 nm、5 nm 步长（81波长），输出1201×1201、角域±6°、接收距离10 m。有限圆盘光源默认每轴5点，光源类型仍由项目选择，默认plane。当前反射分支采用指定频率网格的二次相位卷积傅里叶积分；`paddingFactor` 保留兼容，不再控制反射输出采样。`cfg.heightMap=0`；只计算和展示一组反射结果，膜层自身的相位始终保留。两个单图层校验项目保留小尺寸输出，光谱同样使用5 nm。

<a id="project-config"></a>

## 项目配置 config.txt

配置在光谱、几何和传播计算前加载，**覆盖主脚本默认值**。以下示例应写入项目 config.txt；设置 `useProjectConfig=false` 时才完全采用主脚本参数。完整常用配置已经提供在 Case1~3 中。

常改参数如下（仅示例，完整配置已放入三个 Case）：

```text
# 物理尺度：单位 m。所有图片覆盖完整物面，不自动缩放和平铺。
cfg.pitch = 50e-6
cfg.periods = [16, 16]
cfg.samplesPerPeriod = 64

# 图片文件名可不写扩展名
stack.layerNames = ["Anode", "IJP", "BM", "CF_R", "CF_G", "CF_B"]
stack.layerLevels = [1, 2, 3, 4, 4, 4]
stack.layers.Anode.file = "anode"
stack.layers.Anode.method = "reflectance"
stack.layers.Anode.reflectance = 1
stack.layers.Anode.phaseDeg = 180

# 材料 nk 和厚度；文件非空时优先于 n/k 常数
stack.layers.BM.nkFile = "Materials/BM_nk.csv"
stack.layers.BM.thickness = 1e-6

# 光谱、观察角域与显示
spectrumMode = "D65"
wavelengthStepNm = 5
observation.halfAngleDeg = 2.5
colorOptions.exposureEV = 0
colorOptions.toneStrength = 3000
displayStyle.showAnnotations = true
cfg.heightMap = 0
```

### 统一的图层声明与两种光学设定

不再选择 simple/multilayer。`stack.layerNames` 同时决定读取哪些图片/材料；未声明的图层不读取。名称必须是唯一的有效 MATLAB 字段名，可自定义。`stack.layerLevels` 与列表一一对应，正整数按从下到上非递减排列；省略时按列表顺序逐项叠放。相同层级表示同一物理层中的互斥材料区域，其图片不允许交叠。因此默认 CF_R/G/B 均为层级4，绝不是R膜、G膜、B膜依次上下叠加。

每项可设置 `file`（默认同名图片，可省略扩展名）、`field`（MAT字段名，默认同名）、`invert`（默认false，白色=材料；白色=开口时true）和以下两种 method：

| method | 使用参数 | 物理含义 |
|---|---|---|
| nk | n、k或nkFile，以及thickness(m) | 有限厚度相干膜，包含前后界面及内部多次干涉 |
| reflectance | reflectance(0–1)、phaseDeg(°)，thickness=0 | 不透明有效反射边界，r=sqrt(R)*exp(i*phase) |

直接反射率只有R不足以确定复系数，默认相位180°，可按已知材料指定；R与相位可为常数或逐波长向量。该方式不提供透射系数，不能据此模拟半透明反射片；其下方的光学结构在该位置被遮蔽，其上方nk膜仍按完整膜堆计算。未覆盖区域继续计算其它已声明图层，绝不扩大阳极范围。nk 模式不采用同项的直接反射率参数。

无基底有限膜使用空气作为出射介质：`stack.backgroundMode="medium"`、`stack.incidentN=1`、`stack.exitN=1`。此时未覆盖处无界面反射，入射光可透射；程序当前输出反射量，不把1-R自动报告为吸收。COE项目用 `backgroundMode="absorber"`，未覆盖反射边界的位置为理想匹配吸收终端；两者都可以得到裸底r=0，但物理解释不同。

单图层配置例子（复制到新项目config.txt，放入singlelayer.bmp，并设置匹配图片的cfg采样参数）：

```text
stack.layerNames = ["singlelayer"]
stack.layerLevels = [1]
stack.backgroundMode = "medium"
stack.incidentN = 1
stack.exitN = 1
stack.layers.singlelayer.method = "nk"
stack.layers.singlelayer.n = 1.5
stack.layers.singlelayer.k = 0
stack.layers.singlelayer.thickness = 100e-9
```

改为直接反射率时替换 method，并设置：

```text
stack.layers.singlelayer.method = "reflectance"
stack.layers.singlelayer.thickness = 0
stack.layers.singlelayer.reflectance = 0.36
stack.layers.singlelayer.phaseDeg = 180
```

可直接选择项目 CaseSingleNK 或 CaseSingleReflectance，分别验证上述两种设定；两者均只读取singlelayer.bmp。旧项目需将 stack.materials.* 改成 stack.layers.实际名称.*，删除opticalModel及固定anodeFile/bmFile/cfFiles字段。Case1~3已迁移。

### 光源空间形式（不模拟相机）

屏幕本身不发光，程序计算照明被膜堆反射/衍射后的理想角分布及接收平面辐照度。光源空间形式与光谱独立选择：

| cfg.source.type | 物面入射场 | 参数及定义 |
|---|---|---|
| plane | 理想平面波 | 仍受计算面/光瞳范围限制，保持旧计算行为 |
| point | 有限距离球面波；可扩展为均匀发光圆盘 | distance为法向距离(m)，diameter为发光直径(m)，centerXY为中心投影(m)，beamHalfAngleDeg为辐射强度降至轴向一半的半角 |
| laser | 近轴基模TEM00高斯场，保留波前曲率 | waistRadius为束腰1/e²强度半径(m)，waistDistance为沿传播方向从束腰到物面的有符号距离(m) |

有限圆盘各源点之间非相干，**各自传播后叠加强度，不能将不同源点复振幅直接相加**。samplesPerAxis默认5，使用圆盘内的确定性二维采样，当前5×5候选网格在圆盘内实际保留21点；diameter=0时只用一个点。尺寸非零时至少3×3，当前阶段按用户选择统一使用5×5；计算量随源点数增加。

点源示例（更改同名已有key，勿重复追加）：

```text
cfg.source.type = "point"
cfg.source.distance = 0.1
cfg.source.diameter = 0.2e-3
cfg.source.centerXY = [0, 0]
cfg.source.samplesPerAxis = 5
cfg.source.beamHalfAngleDeg = 30
```

点源角度包络使用轴对称高斯辐射强度分布，beamHalfAngleDeg是半强度半角，不是硬截止角；同时包含几何距离和入射面投影因子。激光并非都属于高斯光束，本程序先支持常见基模；其它横模暂不模拟。

```text
cfg.source.type = "laser"
cfg.source.laserSpecification = "waist"
cfg.source.waistRadius = 100e-6
cfg.source.waistDistance = 0
spectrumMode = "monochromatic"
monochromaticWavelengthNm = 532
```

理想高斯束的远场1/e²发散半角θ=λ/(πw0)，瑞利长度zR=πw0²/λ，束宽w(z)=w0*sqrt(1+(z/zR)²)。束腰与发散角不能自由独立指定。可改laserSpecification="divergence"并填写divergenceHalfAngleDeg，由所选波长反推束腰；此时waistRadius不参与计算。当前限近轴基模，发散半角≤5°。[Newport 高斯束说明](https://www.newport.com/n/gaussian-beam-optics)

所有光源的总入射功率标定到 **cfg.incidentIrradiance × 几何光瞳面积**：incidentIrradiance表示物面平均辐照度，不是激光峰值或LED的总辐射功率。因此改变距离/尺寸会改变波前与照明分布，程序不会凭未知LED功率推算绝对照度的反平方变化。result.incidentIrradianceMap保存代表波长的物面照明，result.illumination保存各波长模式数、归一化因子、束宽和发散角。

输出传播仍为夫琅禾费角分布/接收平面，距离不足时有远场提示。膜堆仍采用法向局部复系数；源导致局部入射角超过5°时提示近似限制，超出物面采样带宽则报错。不包括相机、斜入射偏振修正或电磁侧向耦合。

### 光谱：D65、自定义与单色

- spectrumMode="D65"：使用现有CIE D65；默认380–780 nm、5 nm步长。
- spectrumMode="custom"：sourceSpectrumFile指向项目中的两列CSV（无表头），依次为波长nm、非负相对光谱功率密度。波长严格递增，覆盖wavelengthRangeNm，不外推；线性插值到计算网格，使用梯形积分。不能将RGB值或色温填入此文件。
- spectrumMode="monochromatic"：按monochromaticWavelengthNm计算单条谱线，总谱线功率归一化为1，不把它当作1 nm宽的连续谱。
- RGB保留旧三波长代理预览；前三种模式均使用CIE三刺激值函数合成色彩图，保持原始XYZ及物理指标。当前色度范围为380–780 nm。

```text
spectrumMode = "custom"
sourceSpectrumFile = "my_source.csv"
wavelengthRangeNm = [380, 780]
wavelengthStepNm = 5
```

文件格式例子（以下仅为格式示意，数值应替换为自己的实测光谱）：

```csv
380,0.1
450,1.0
550,0.6
650,0.3
780,0.05
```

D65沿用*_d65_color.png；自定义/单色输出*_spectral_color.png，*_spectrum.csv记录实际IncidentSPD。源SPD只定义相对光谱形状，功率标定仍由物面平均辐照度决定；辐射度和明视觉分母分别按实际入射光谱归一化。

### IJP：统一用全白图片输入

默认顺序是 **Anode → IJP → BM → [CF_R,CF_G,CF_B]**。IJP全平面连续覆盖，上下表面平坦，n=1.5、k=0、厚10 μm。Case1~3已提供与其它图层严格同尺寸、同原点的全白IJP.png/.bmp/.jpg，不再需要特殊的无图片输入流程。白色表示“此处有膜”，透明性由k=0定义；透明膜仍有界面反射和传播干涉，绝不能以删除复反射相位的方式模拟透明。

```text
stack.layers.IJP.file = "IJP"
stack.layers.IJP.method = "nk"
stack.layers.IJP.n = 1.5
stack.layers.IJP.k = 0
stack.layers.IJP.thickness = 10e-6
```

IJP下方仅Anode图样范围内为反射底，其余位置沿用理想匹配吸收终端；覆盖IJP/BM/CF后的界面反射仍计算。工程接口保留coverage="full"便于算法验证，默认项目统一读取全白图片。

厚膜干涉对光谱步长更敏感：仅10 μm、n=1.5的IJP在380 nm附近，往返相位变化周期约4.8 nm，5 nm可能不能分辨局部干涉条纹。当前常规方案比较按用户选择暂用统一5 nm；对很接近的候选排序或精细光谱形状，再用更细步长复核。主脚本保留厚膜采样提醒；完整膜堆的最严要求可能比单独IJP更高，建议步长不是精度保证。

### 采样收敛基准

选定项目后，将主脚本 `runSamplingBenchmark=true`，完整运行。默认关闭；开启后复用项目配置、用户图片与统一膜堆，完成基准后返回，不执行常规绘图/导出。实现位于 `src/private/benchmarkCoeSampling.m`，不增加公共接口。基准要求图片输入、零额外高度、空外部复系数及无限光瞳；不调用archive，不改写Input。

这是独立的固定条件实验：D65 380–780 nm、平面照明、100 m夫琅禾费接收距离；有限源实验另设距离0.1 m、直径0.2 mm的均匀非相干发光圆盘。参数与项目日常配置分开记录。首个基准使用Case1：50 μm Pitch、16×16单元、每Pitch 64点。其它尺度应先核对Nyquist角域、ROI完整性及计算资源。

| 扫描 | 设置与参照 | 判断范围 |
|---|---|---|
| 光谱 | 5/2/1/0.5/0.25/0.125 nm，相对0.0625 nm；另试半格偏移并保留波段端点 | 区域R曲线、宽带总反射率、固定(m,n)阶次效率及原始XYZ/xy |
| 单独IJP | 同一全白图片，空气/10 μm、n=1.5、k=0/空气 | 与独立Airy公式核对，区分单膜与完整膜堆影响 |
| 输出角网格 | ±1°：201/401/801/1601；带宽允许时另测±6°：401/801/1201/1601 | 380/550/780 nm；最少每像素每轴3点积分 |
| 像素内积分 | 固定±1°、401输出，最少1/3/5/7点 | 同时记录自动提高后的实际点数 |
| 重复单元 | 每轴8/16/32/64个Pitch，固定物面步距、±1°/801输出 | 所有图层逐像素验证完整2×Pitch周期后，仅在基准内存中显式平铺；非周期输入跳过 |
| 物面栅格 | 原图共同原点整数抽取至每Pitch 16/32/64点 | 参照原生64点；真实的更细几何收敛需要更高分辨率输入 |
| 硬像素积分 | 每Pitch 64/128点，零阶保持原像素支撑 | 只检验同一栅格图的数值积分，不能恢复亚像素边界 |
| 有限源 | 550 nm、8×8 Pitch、±1°/201输出，每轴3/5/9/17/25点 | 实际圆盘内点数另列；非周期输入使用原全物面；各点分别传播后叠加强度 |

空间扫描比较固定0.8°圆锥、理论(1,0)阶次0.08°/0.10°圆域、(2,0)及(1,1)阶次0.08°圆域，直接按生产输出bin中心判断圆边界，包含ROI离散误差，不跨组追踪S编号。另按bin重叠面积守恒重分箱到共同±1°的200×200网格，报告角分布L1绝对差；共同窗口边界固定，避免输出端点bin宽度变化混入窗口损失。峰值与网格FWHM仅作诊断，周期数增大时峰会变高变窄，不能要求峰值收敛。

初始容差为 `绝对差 ≤ 1e-4 + 1%×参照量`（反射率0–1），分别检查积分量与共同角分布。局部光谱曲线另看最大绝对误差，不能用宽带平均替代。最细样本自比较用 `IsReference` 标记，其通过不代表独立精度证明。重复数改变了物面的真实物理尺寸，差异不能全归为数值误差。

输出在 `Output/项目名称/SamplingConvergence/`：`spectral.csv`、`ijp_spectral.csv`、`spatial.csv`、`convergence.png`；`spectral_reference.mat`保留区域/阶次曲线和XYZ；`sampling_convergence.mat`保留原始regionMap、区域/材料表、参数、共同角图、项目配置文本及源码快照。`spatial_partial.csv`随扫描更新。派生物面由原图与Family/Value重建，不写入用户图片目录。

光谱扫描通过区域系数及空间傅里叶系数避免数千波长的全物面传播数组。严格周期时阶次功率代表无限阵列；非周期时仅为指定空间频率诊断。它不替代有限物面固定ROI或完整D65角图的联合光谱收敛验收；空间扫描先用单色隔离因素。Airy解析与离散Parseval核对也不替代接收面积分收敛。

日常不同方案的初步比较统一采用5 nm，保持光谱网格、光源和其它采样条件相同；不要求所有方案改用0.25 nm。局部光谱曲线的最大重建误差不能直接当成方案差值或排序误差。对差异很小的候选，再以更细步长抽查结论是否稳定。

“有限源采样”指照明光源发光面的积分，不是版图重复数。例如每轴17点是在发光圆盘的包围方形内建立17×17候选点，本次圆盘内实际225点；每点分别照明完整版图，再叠加强度。版图重复数由cfg.periods独立控制，本次有限源扫描为[8,8]，常规Case1为[16,16]。

`angle_grid_380nm.png`、`angle_grid_550nm.png`、`angle_grid_780nm.png`提供直接图片对照：左为±6°/801全图，中为同一结果裁切显示±1°，右为真正按±1°/801重新计算。中栏只有约133个原始角度样本横跨当前视窗，右栏有801个；两者间距分别为0.015°与0.0025°。三栏使用共同能量密度dB标尺；均为平面波、16×16 Pitch、单色空间采样对照，不是D65彩色合成。原始数组和显示参考保存在`angle_grid_comparison.mat`。

2026-11-13首次Case1实跑完成68组空间扫描。下表的光谱误差是相对0.0625 nm曲线、使用线性重建的最大R绝对差，取两个采样起点的较差值；单位为**百分点**，不是相对百分数：

| 步长 | 单独10 μm IJP | 完整Case1膜堆最敏感区域 |
|---|---:|---:|
| 5 nm | 14.67 | 30.05 |
| 1 nm | 1.65 | 5.06 |
| 0.5 nm | 0.43 | 1.38 |
| 0.25 nm | 0.11 | 0.37 |

5 nm时，完整膜堆的D65总辐射反射率约9.7200%（半格偏移后9.6076%），参照为9.6683%；总量看起来接近，仍不能代表干涉条纹已分辨。1 nm的宽带总量误差已经很小，但精细光谱对照应从0.25 nm及以下检查；这不保证尖锐自定义SPD或固定角ROI也已收敛。单独IJP生产模型与独立Airy公式的最大差为1.11e-16。

- ±1°/801相对1601的共同角图L1差占参照窗口能量约0.27–1.01%；±6°/801相对1601为9.04–25.70%，而所选积分指标最大差仅约0.0017个百分点。默认宽视图仍不宜用于判断细峰形状已经收敛。
- 固定±1°/401时，像素内最少3点相对7点的共同角图差约0.044–0.154%；该组条件下默认3点满足初始混合容差。
- 32→64个Pitch的最大积分变化约0.0122–0.0188个百分点，780 nm仍有阶次ROI未达初始容差；峰形随物面尺寸变化更明显，不能直接宣布有限周期已充分。
- 每Pitch 32点相对原生64点的共同角图差约1.99–3.21%；64→128硬像素保持细化的差不超过约0.108%。后者说明现有栅格的积分较稳定，不能证明原生几何分辨率已经充分。
- 有限圆盘3×3实际9源点，相对25×25的489源点，共同角图差约6.40%；17×17实际225源点降至约0.070%。5×5、9×9误差不单调，来自圆盘边界的离散采样；不能仅凭一对粗网格宣布收敛。该结论限于本次8×8 Pitch、550 nm及指定光源尺寸/距离。

68组的最大离散Parseval误差为1.31e-14；现有37项COE回归通过。以上是示例材料和图样的数值基准，不是产品实测精度或所有采样维度联合收敛的保证。

### 每轴需要多少Pitch：专项基准

在主脚本中设置 `runSamplingBenchmark=true`、`samplingBenchmarkMode="periods"`。原五类扫描使用mode="all"。专项读取同一项目图片及统一膜堆，逐像素确认完整2×Pitch晶胞后，扫描每轴8、16、24、32、48、64、96、128、192、256个Pitch；结果写入 `Output/项目名称/PitchConvergence/`。这里的128指128×128个Pitch，沿每轴64个完整2×Pitch晶胞。

本项固定平面波、D65 380–780 nm/5 nm、±6°/1201输出；每Pitch保持原生64点。以固定窗口能量及理论阶次的固定方形ROI效率判断收敛：零级、(1,0)、(0,1)、(1,1)、(2,0)、(0,2)、(2,2)、(4,0)、(4,4)、(8,0)、(12,0)，每个分别使用0.08°和0.10°半宽。两种半宽均不覆盖相邻阶次。此处明确使用**方形ROI**，不能与上一轮圆形ROI表逐行混比，也不把峰高或窄峰逐点形状用作尺寸收敛条件。

专项通过晶胞复场的非周期自相关计算有限阵列强度，晶胞间距的相关项权重为`R-abs(t)`（R为每轴完整晶胞数），对每个角度矩形解析积分。它与生产传播采用相同的离散物面点求积，不增加像素孔径、不删膜堆相位；避免展开8192×8192以上的大物面，也避免窄峰细分网格造成的巨大数组。另用现有runCOEReflection入口在380/550/780 nm、8/16 Pitch做6组输出图交叉检查，并对每个波长/重复数独立核对整个Nyquist方窗的Parseval。该加速只在基准中使用，正式输入及传播流程保持原有规则。

参照为同一晶胞DFT给出的**无限周期阵列阶次效率**，比仅比较相邻有限尺寸更严格。沿用容差`绝对差 ≤ 1e-4 + 1%×参照量`，其中1e-4是0.01个百分点；分别报告D65辐射度/明视觉和81个逐波长是否通过。推荐点须由至少一个更大样本继续通过确认，不能用最末样本自比较为零宣布收敛。

2026-11-13 Case1结果（误差为所测ROI/窗口中的最大绝对差，单位**百分点**）：

| 每轴Pitch | 物面边长 | D65辐射度最大差 | D65明视觉最大差 | D65均通过 | 81波长均通过 |
|---|---|---:|---:|---|---|
| 16 | 0.8 mm | 0.22948 | 0.21517 | 否 | 否 |
| 32 | 1.6 mm | 0.11879 | 0.11026 | 否 | 否 |
| 64 | 3.2 mm | 0.05957 | 0.05574 | 否 | 否 |
| 96 | 4.8 mm | 0.03980 | 0.03723 | 否 | 否 |
| 128 | 6.4 mm | 0.02989 | 0.02796 | 是 | 否 |
| 192 | 9.6 mm | 0.01996 | 0.01866 | 是 | 是 |
| 256 | 12.8 mm | 0.01498 | 0.01400 | 是 | 是 |

在已测档位中，**常规D65/5 nm方案比较可先用每轴128 Pitch**；若还要求81个单波长的所测ROI全部满足同一容差，使用**每轴192 Pitch**，并由256 Pitch确认。最敏感项是0.08°半宽的零级积分，有限尺寸旁瓣泄漏随阵列增大缓慢减少；只看总反射率或相邻两组变化会较早误判收敛。128 Pitch的该项D65反射率为2.22797%，无限阵列为2.25786%。

这是当前平面波、图样、材料及所列ROI的结论，不是所有光源/观察条件的通用Pitch数。有限距离点源的入射曲率和照明范围会随物面尺寸改变；选择5×5光源采样不消除这种物理变化。正式改用128/192 Pitch时，应准备对应的8192×8192/12288×12288完整同坐标用户图片，且核对远场条件与资源；程序未自动更改现有16-Pitch输入配置。

`pitch_summary.csv`保留逐波长和宽带的分开通过标志及推荐；`pitch_convergence.mat`保留各波长/尺寸/ROI原始效率、无限阵列参照、完整配置与源码快照；`production_validation.csv`保存交叉验证；`pitch_convergence.png`给出误差曲线，`pitch_images_550nm.png`在共同dB标尺下展示不同Pitch数的单色角图。当前交叉验证最大角图L1相对差约0.095%（生产代码采用有限midpoint积分），精确周期积分的最大Parseval差约7.5e-15。

### 显示范围与默认细腻程度

照片未给角度标尺，不能从照片像素直接判断±6°是否覆盖足够级次。当前采用±6°全局角域、outputSize=1201，角间距0.01°。这是真实输出采样，像素内仍自动细分积分。1201相比旧801需要约2.25倍的光谱数组存储，计算时间也可能增加。

主脚本会报告角间距、最短波长的λ/D特征尺度及其包含的显示间隔。50 μm Pitch、16×16单元、800 μm物面在380 nm时λ/D约0.0272°，1201网格约2.7个间隔；增大Pitch数后峰进一步变窄，因此积分效率通过不意味着图像分辨了真实峰宽。具体范围仍应结合CapturedFraction和光源参数判断。

若主要看中心细节，可使用真实局部计算（不是只改colorView裁切）：

```text
cfg.angleHalfRangeDeg = 3
cfg.outputSize = 1201
```

角间距此时0.005°，仍可覆盖2.5°观察锥。若改为更小窗口，需同步减小observation.halfAngleDeg、cfg.metricAngleDeg及其它越界ROI。扩大到±10°时，为保持0.015°间距可用outputSize=1335，并先核对物面Nyquist角域及近轴限制。仅改变colorView.halfAngleDeg会裁切原图，不会增加真实细节。先用单色或RGB代理检查采样，再用D65/自定义全光谱做最终比较，避免一次性增加空间、源点和波长采样导致内存/时间陡增。

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

当前使用正方形采样像素。`periods` 采用偶数与现有配置保持一致。若导出的图样每个 50 μm 单元有 100 点，则将 `samplesPerPeriod` 改为 100，并按上述公式核对完整矩阵尺寸。需要保留小于原采样间距的窄边/间隙时，应先提高真实版图的栅格采样；当前反射传播由输出角域与 `outputSize` 控制接收像素，`receiverBinSamples` 控制像素内积分，`paddingFactor` 不再控制该分支，三者均不能修复物面欠采样。

### 图片模式下比较 16、32、64 个 Pitch 单元

保持 `cfg.pitch=50e-6`、`cfg.samplesPerPeriod=64`，则每像素始终为 0.78125 μm。六张图 anode、IJP、BM、CF_R/G/B 必须同尺寸、同原点、同样复制；优先使用无损二值 PNG，不能用图片缩放代替增加重复单元。

| 项目 | cfg.periods | 每张图片尺寸（行×列） | 完整物面尺寸 |
|---|---|---|---|
| Period16 | [16,16] | 1024×1024 | 800×800 μm |
| Period32 | [32,32] | 2048×2048 | 1600×1600 μm |
| Period64 | [64,64] | 4096×4096 | 3200×3200 μm |

这里的计数是 Pitch 单元数，未必是完整光学晶胞数。Diamond 的轴向完整晶胞为 2*Pitch，因此三组沿每轴分别含 8、16、32 个完整晶胞。应以包含各层材料、厚度及高度相位的复反射分布的真实平移周期判断级次。

准备方法：先确认一个完整晶胞（Diamond 示例为128×128像素）在各图层中均首尾周期衔接，然后分别铺8×8、16×16、32×32次。若已有严格周期且边界正确的1024×1024全物面图，也可原样保留16组，并将整张图分别平铺2×2、4×4次生成32和64组。裁取晶胞后不得重复包含首尾同一个采样端点；各层必须使用相同裁切位置。完整16组中若有缺陷或非周期边缘，平铺会制造更大的缺陷超晶胞，不等同于目标无限光栅。

各项目复制相同的 Materials 和 config.txt，修改 periods 与对应图样路径，主脚本用 projectName 切换。对边界衔接正确的无损 Case1 图，可用以下显式离线操作准备全物面图片（在模块目录执行；不要覆盖现有项目）：

```matlab
layers = ["anode", "IJP", "BM", "CF_R", "CF_G", "CF_B"];
for count = [16, 32, 64]
    target = fullfile("Input", "Period" + count);
    assert(~isfolder(target), '目标项目已存在，请换名，避免覆盖输入。');
    mkdir(target);
    copyfile(fullfile("Input", "Case1", "Materials"), fullfile(target, "Materials"));
    for layer = layers
        image = imread(fullfile("Input", "Case1", layer + ".png"));
        assert(size(image,1)==1024 && size(image,2)==1024);
        imwrite(repmat(image, [count/16, count/16, 1]), fullfile(target, layer + ".png"));
    end
    text = fileread(fullfile("Input", "Case1", "config.txt"));
    text = strrep(text, 'cfg.periods = [16, 16]', ...
        sprintf('cfg.periods = [%d, %d]', count, count));
    writelines(string(text), fullfile(target, "config.txt"), 'Encoding', 'UTF-8');
end
```

先用三波长代理（`spectrumMode="RGB"`）检查周期收敛，再用 D65 光谱检查正式结果，避免宽带色散模糊掩盖单波长窄峰。三组使用相同波长、材料、厚度、入射辐照度、距离及输出角域。可在各组 config 中统一改为以下局部窗口配置；已有同名 key 必须替换，不能重复追加：

```text
spectrumMode = "RGB"
cfg.propagationDistance = 100
cfg.incidentIrradiance = 1
cfg.angleHalfRangeDeg = 1
cfg.outputSize = 1201
cfg.receiverBinSamples = 3
cfg.metricAngleDeg = 0.8
observation.halfAngleDeg = 0.8
colorView.halfAngleDeg = 0.8
spotOptions.integrationHalfAngleDeg = 0.08
spotOptions.zeroOrderExclusionDeg = 0.10
spotOptions.minSeparationDeg = 0.18
```

`pupilRadius` 保持默认 Inf，heightMap 保持0；若使用有限光瞳，需使有效照明宽度随周期数一起增大，否则增加未照明区域不改变结果。100 m 在64组、380 nm下大于常用远场尺度2*D²/λ≈54 m；这是自由空间夫琅禾费近似的距离设置，不表示实际装置必须放置100 m接收面。若实际接收距离较短，应评估透镜傅里叶面或近场模型，不能靠增加z人为解释原测量。

上述±1°/1201网格的角间距为0.001667°，64组在380 nm的峰中心至第一零点约0.0068°，有约4个间隔；比默认±6°/401的0.03°更适于分辨窄峰。该局部窗口不覆盖全部级次，出现窗口捕获不足98%的提示是预期现象。若保留±6°且要求同样角度分辨率，需要7201点/轴，计算和存储会明显增加。自动像素内细采样可以改善积分，不能使粗像素图显示出窄峰的真实宽度。

收敛比较以逐波长的理论级次方向为固定中心，使用三组完全相同的积分半角；应覆盖64组和16组的峰及足够侧瓣，并保持相邻级次ROI不交叠。将积分半角0.08°与0.10°作敏感性比较；若边缘级次ROI越界则扩大输出窗口或不比较该级次。S1/S2是自动亮斑排序号，不是固定(m,n)编号，不能直接跨组按编号匹配。比较同一角位置的绝对反射率，并观察32→64的变化是否小于所需容差，例如相对1%（弱级次同时设绝对容差）。峰应随尺寸增大而变窄，峰高不必收敛；各组独立自动增强的彩色图不能用来判断效率收敛。先在固定周期数下增加samplesPerPeriod并重新栅格化图片，验证物面采样收敛；再在固定物面下增加outputSize、receiverBinSamples验证输出积分收敛，才能把剩余变化归因于有限周期数。

### 图片或 MAT 图样输入

默认 `stack.geometry.source="images"`，支持 PNG、BMP、JPG/JPEG、TIF/TIFF。推荐无透明通道/调色板的灰度PNG/BMP；JPEG有损压缩可能影响边界。所有图层严格同尺寸、同原点，图中不得包含标尺/文字。灰度归一化后按 `stack.geometry.threshold`（默认0.5）二值化；RGB图片取通道平均，不能用一张合并的彩色像素图代替各材料的独立掩膜。

每个列表项的 `stack.layers.名称.file` 指定图片，省略时使用名称自动匹配扩展名。同名多格式时需明确扩展名。默认 BM 图片白色=实心材料；若实际图片白色=开口，设置 `stack.layers.BM.invert=true`。

MAT模式配置如下（config.txt语法）：

```text
stack.geometry.source = "mat"
stack.geometry.matFile = "layers.mat"
stack.layers.Anode.field = "anodeMask"
stack.layers.BM.field = "bmMask"
stack.layers.CF_R.field = "cfRMask"
stack.layers.CF_G.field = "cfGMask"
stack.layers.CF_B.field = "cfBMask"
```

只需提供已声明列表项的顶层字段，值为完整Ny×Nx的logical或严格0/1矩阵。不同项目可以有不同名称与层数，未声明的字段不参与计算。若提供bmOpenMask，设置BM.field="bmOpenMask"和BM.invert=true。程序不再自动外扩或生成缺失图层，BM/CF轮廓以输入图样为准；同层材料不交叠，跨层图样可以交叠。

### 材料 nk 数据与厚度

Case1、Case2、Case3 的 `Materials/` 目录也提供 `CF_R_nk.csv`、`CF_G_nk.csv`、`CF_B_nk.csv` 样例，均覆盖 380–780 nm。R 使用 n=1.65、k=0.03，G 使用 n=1.60、k=0.04，B 使用 n=1.55、k=0.05；均沿用程序默认占位常数，**不是实测数据，也不表示真实滤色片的光谱特性**。文件格式与 BM 相同，可直接替换为实测数据；在项目 `config.txt` 中按下方对应路径启用，厚度仍单独设置。

`Input/Case1/Materials/BM_nk.csv`、`Input/Case2/Materials/BM_nk.csv` 和 `Input/Case3/Materials/BM_nk.csv` 已提供可读取的样例，覆盖 380–780 nm，沿用默认占位常数 n=1.80、k=0.35，**不是实测数据**。后续可直接按同样格式替换。样例未默认启用；使用时将对应项目 `config.txt` 中的 `stack.layers.BM.nkFile` 改为 `"Materials/BM_nk.csv"`。默认启用项目配置时，配置中的值会覆盖主脚本中的设定；该路径相对于当前项目文件夹。

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
stack.layers.BM.nkFile = "Materials/BM_nk.csv";
stack.layers.CF_R.nkFile  = "Materials/CF_R_nk.csv";
stack.layers.CF_G.nkFile  = "Materials/CF_G_nk.csv";
stack.layers.CF_B.nkFile  = "Materials/CF_B_nk.csv";

% 以下厚度数值只是填写示例，请换成各层实际厚度；单位 m。
stack.layers.BM.thickness = 1000e-9;
stack.layers.CF_R.thickness  = 1100e-9;
stack.layers.CF_G.thickness  = 1000e-9;
stack.layers.CF_B.thickness  = 900e-9;
```

该层 `method="nk"` 时，`nkFile` 非空会优先于该材料的 `n`、`k` 常数/向量；厚度仍由 `thickness` 单独指定，不写在 CSV 第四列。`nkFile=""` 则使用脚本中的 n/k，可为常数或与计算波长一一对应的向量。

默认Anode采用method="reflectance"、R=1、phaseDeg=180、thickness=0；也可以将其改为nk有限膜。各区域是否存在只由其图片决定。

### 检查和输出

自动检查名称、顺序、图片尺寸、材料范围以及同层不交叠；程序不替用户修改BM/CF图片或补齐覆盖区域。COE图片应继续符合阳极/BM开口/CF覆盖的物理约定。

运行前可设 `caseName="actual_coe_01"` 区分输出；同名文件会覆盖。优先核对输出几何图、`*_layer_patterns.mat`、`*_stack_materials.csv`，确认最终掩膜、插值后的 n/k 和厚度。然后查看 D65 图、累计角能量与绝对反射率 CSV。显示增强不改变原始 XYZ 和物理指标。

额外高度保持 0 nm 时只计算一次反射并跳过额外高度相位；膜层自身的厚度干涉相位始终保留。当前仍是法向局部相干薄膜 + 标量近轴远场模型，忽略侧壁和横向电磁耦合。

<a id="model"></a>

## 光学模型与计算流程

主脚本按“初始化 → 读取项目配置 → 读取图样/材料 → 反射计算（非零高度才清零并再次计算） → 积分评价 → 展示/保存”顺序组织。具体算法在 private 中，用户入口保持完整、可从上往下阅读。

### 阳极、BM 与 CF

- 平面像素采用 Diamond 排布；用户提供 Pitch，例如 50 μm。
- 从下到上为阳极、IJP、BM、RGB CF。各层是图案化区域，并非整面均匀覆盖。
- **阳极反射只存在于输入阳极图样范围内**，其裸底反射率默认 100%。
- 默认示例的BM开口外扩阳极5 μm；正式计算直接采用输入BM图样，开口外为实心BM。
- **BM 开口中超出阳极图样的外扩环带属于不透明吸光底层，不能按 100% 阳极处理。**
- R/G/B CF 覆盖各自对应的 BM 开口，其区域大于 BM 开口；各色 CF 区域之间不得交叠。
- 用户可输入阳极、BM 和 RGB CF 图样，并提供 BM / R CF / G CF / B CF 的 n、k 与厚度。

### 局部膜堆

| 区域 | 从入射侧向下的局部结构 |
|---|---|
| 阳极内、BM 开口内、有 CF | 空气 / 对应 CF / IJP / 阳极反射底 |
| 阳极外、BM 外扩开口内、有 CF | 空气 / 对应 CF / IJP / 吸光底层 |
| BM 实体上方有 CF | 空气 / 对应 CF / BM / IJP / 该位置的底层 |
| BM 实体上方没有 CF | 空气 / BM / IJP / 该位置的底层 |

底层由阳极掩膜判定，不能仅依据 BM 开口判定。阳极默认裸底 r=-1（R=100%、相位180°），反射率和相位均可调。阳极外的理想吸收终端在直接接触入射介质时裸底 r=0，即阻抗匹配；进入终端的净能量被吸收，不作为背后透射光再次返回。加上 CF/BM 后，按同一个终端阻抗重新计算膜层/终端界面的反射，所以局部最终系数不一定为零。这是理想边界，不是已经测得的吸光材料 n/k。

输入图样必须有一致的物理尺度、坐标原点和矩阵尺寸。程序只按声明读取并计算各区域；同一层级的区域必须互斥。实际COE图片仍应满足BM开口包含阳极、CF覆盖对应开口等约定。

各位置先按局部平行平面膜堆计算随波长变化的复反射系数，再进入指定频率网格的傅里叶积分传播 / 角域 / CIE 后处理。CF 局部厚度固定，位于 BM 实体上方时叠加在 BM 顶面，在开口内则从裸底平面起算。以各物理层级最大nk厚度之和为共同空气参考面（默认COE即 `IJP厚度+BM厚度+最大CF厚度`），补足各柱体顶面到该平面的往返相位。未包含侧壁或邻近柱体的矢量耦合，不能当作 RCWA/FDTD。

膜厚传播与干涉相位属于膜层模型，不等同于可选的额外高度图；`heightMap` 默认仍为 0。几何采用 regionMap，局部光谱系数采用实际存在的结构数×波长数的查找表，避免构建庞大的逐像素光谱系数三维数组。



该模型为法向局部平行平面相干膜堆 + 标量近轴远场，已支持有限圆盘源的非相干强度积分；仍忽略侧壁、邻近柱体矢量耦合、斜入射膜堆修正及偏振，不是RCWA/FDTD或真实全半球散射。常数CF n/k不保证实际RGB通带，需输入真实色散吸收数据。

### 光谱与色度积分

令 S_i 为 D65 相对光谱、q_i 为梯形积分宽度（nm）。各波长的输出数组已相对该波长入射能量归一化。

- 辐射度：`w_i = S_i*q_i / sum(S*q)`。
- 明视觉：`wY_i = S_i*q_i*ybar_i / sum(S*q*ybar)`。
- XYZ：`WXYZ_i = S_i*q_i*[xbar_i,ybar_i,zbar_i] / sum(S*q*ybar)`，三个分量共用一个归一化常数。

因此辐射度和明视觉绝对反射率分别以入射辐射能与入射明视觉光通量为分母。二者不能混为一项，也不能先对每个波长的“反射能量内部比例”求平均来替代光谱能量积分。

所有光谱量只覆盖 380–780 nm，未把波段外的入射能量计入分母。D65 在该范围、5 nm 步长下得到白点 xy 约 (0.31272, 0.32904)，与 D65 标准色度接近；这是数据及积分的检查，不代表波长步长已经对所有尖锐衍射结构收敛。

### 两种高度状态与兼容入口


`heightMap` 为标量0或全零矩阵时自动采用单状态模式，无需另设开关：不计算额外高度相位，不重复传播，也不重复进行累计角能量、光斑状态积分或D65色度合成。远场角度图与接收面图单行排列，D65图横向展示“全范围/局部放大”，指标曲线只绘一条，输入几何图不显示全零高度图。膜层nk、厚度干涉和共同参考面相位仍完整计算。

非零高度（包括均匀非零高度）保留原来的两状态比较。保存的MAT增加heightPhaseEnabled，result.propagationRuns为1或2。baseline、baselineKernel/baselineMetrics、withoutHeight/withHeight等字段继续保留作兼容，单状态时复用同一份数据，不代表执行了第二次计算。

“不叠加高度相位”仅将额外 heightMap 清零，“叠加高度相位”采用输入的 heightMap。两者仍保留膜厚干涉和共同参考面的相位；它们是同一输入结构的两种计算状态。`baseline` 字段名为兼容历史接口保留。

额外高度相位 `phi=4*pi*deltaIndex*h/lambda`，h单位m。均匀高度只增加整体相位，非均匀高度才重新分配远场能量；`deltaIndex` 是额外光程层的有效折射率差，不能直接当作某一种膜的n。不要把膜厚相位重复加入 heightMap。

高度条纹自动生成器已存档，主脚本不再提供useHeightExample/heightExample参数。若需评估非零高度，请由用户准备heightMap矩阵（m），通过inputFile指定的MAT文件加载。旧config中的useHeightExample或heightExample.*应删除；历史0/687.5 nm条纹仅保留在附录和归档源码中。

`inputFile` 是旧单层/额外高度/外部复系数 MAT 入口，与多层图样的 `stack.geometry.matFile` 不同。可包含 mask、heightMap、complexCoefficient：mask为二维0–1区域图，高度为标量或完整二维矩阵，复系数可为标量、二维或Ny×Nx×波长数，模不大于1。主脚本请通过stack.layerNames和各层图片输入几何，不使用旧mask。外部复系数非空时优先于多层模型；若其已包含额外高度相位，应设置 heightMap=0。

主脚本没有单层/多层分支，统一按声明构建局部膜堆。底层兼容公共接口仍可接受历史cfg，但不再作为用户主流程或配置选项。`spectrumMode="RGB"` 仍表示三波长代理预览，与层数无关。

<a id="metrics"></a>

## 指标与结果解释

### 传播距离与接收面强度

当前反射核心与 `topic_aperture_fraunhofer` 使用相同的二次相位卷积傅里叶积分算法（分两轴执行以节省内存），在指定频率网格计算F；物面端点不重复，积分面积元为dx*dy。传播按 `U=exp(i*2*pi*z/lambda)/(i*lambda*z)*F(X/(lambda*z),Y/(lambda*z))`，本次按需求省略观察面二次相位，它不影响强度。仍是夫琅禾费、法向、标量近轴模型，不会因加入z自动成为Fresnel近场求解器。

`cfg.propagationDistance` 为z（m），`cfg.incidentIrradiance` 为每个波长独立标定的入射辐照度（W/m²，默认1），不是D65的每nm绝对光谱辐照度。`result.receiverCoordinate=z*result.coordinate` 为X/Y坐标（m），`result.spectralIrradiance(:,:,k)` 为每个接收像素内平均辐照度（W/m²）；主脚本输出 `*_receiver_irradiance.png`，MAT保留坐标、强度、像素面积、入射面积与入射功率。各波长以同一参考入射值分别计算，不能直接把这些参考强度数组相加当作绝对D65辐照度；绝对宽带输出还需实际入射光谱及积分权重。

长度坐标图与角度坐标的单波长图复用同一组归一化dB数组：参考为两种高度状态、全部计算波长的共同最大值，下限为displayFloorDb（默认-60 dB），上限0 dB，灰度色表和刻度相同。等面积接收像素满足E/Eref=I/Iref，因此这里只变更坐标标签；MAT中的原始spectralIrradiance仍为W/m²，没有被显示归一化覆盖。displayData记录referenceEnergy、等价referenceIrradiance、实际预览波长及两种状态的dB数组。RGB/D65彩色图继续使用其色度显示映射，不转成dB。

每个像素用midpoint细采样积分，`receiverBinSamples` 为每轴最少奇数采样数；程序按物面总宽度和波长自动增加，保证每λ/L至少4点，可进一步提高该参数验证精度。实际细分数记录在 `result.receiverSamplesPerBin`。`spectralEnergy=像素平均辐照度*接收像素面积/入射功率`（入射为零时按单位入射参考计算归一化能量），不会将窗口捕获能量强行归一化。扩大物面可能提高内部积分计算量，即使outputSize保持不变。

固定角度网格、同一物面与入射辐照度时，z加倍后接收面X/Y加倍、辐照度变为1/4，而像素积分能量及反射率不变。固定物理接收尺寸比较不同z时，必须相应缩小angleHalfRangeDeg。程序检查接收像素边界不超过物面采样的Nyquist频带，并在z小于常用远场尺度2*D²/λ_min时提示；该尺度只是近似判据，应根据允许的相位误差判断是否足够远。

`ParsevalAbsoluteError` 仍保留为独立未补零FFT的离散能量核对，只验证离散Parseval关系，不表示指定接收面上的积分误差。接收面积分精度应通过细化receiverBinSamples、输出网格和物面采样分别检查。

人眼距COE 20–50 cm观察反射点光源光晕，与同距离放置裸接收屏的物理量不同：眼睛还通过瞳孔和眼内透镜成像。当前角域结果可用于分析理想平面波照明下的衍射方向和能量分配，但不是视网膜图像。对于距COE 5–50 cm的实际点光源，定量模拟还需球面入射波、实际照明/空间相干范围、COE到眼睛的传播、瞳孔及对焦；自由空间段可采用Fresnel或角谱法，眼内再进行成像。不能仅凭观察距离将整个系统判定为自由空间远场，也不能把propagationDistance直接当作人眼的等效焦距。

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

本模块逐波长的评价公式为 `100*(1 - E_inside_cone/E_total_reflected)`，单位 %。`cfg.metricAngleDeg` 始终表示**半角**，默认 2.5°，圆锥中心是法向入射下的镜面反射方向。分母采用物面反射场的积分能量，并以独立FFT核对Parseval关系，而不是显示窗口内能量，因此指标包含窗口外能量。它是标量近轴模型的估算，也未做标准仪器所需的光度加权及校正。若没有反射能量，该指标为 NaN（未定义）。

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


单状态输出不再重复写入基准/高度对比文件：使用 `*_reflection_metrics.csv`、`*_rgb.png`、`*_srgb_linear.png`；observation_metrics和diffraction_spots仅含一个状态，cumulative_energy的列为RadiusDeg/RadiantPercent/PhotopicPercent。非零高度保留原来的baseline_metrics/comparison_metrics和without/with图片、累计曲线列。两种模式均保留原始MAT、光谱、材料/区域表与公共图像后缀。同名项目目录中旧运行的文件不会自动删除，请以本次配置和MAT中的模式标志判断对应输出。

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
| `_layer_patterns.mat` | layerMasks结构、layerNames和layerLevels；复用为MAT输入时先将layerMasks各字段另存为顶层字段 |
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

### 私有函数审计与生成工具存档

当前 `src/private/` 保留28个函数文件；按主脚本及两个公共接口的静态调用关系检查，没有完全不可达的独立函数。但“可达”包含条件分支，不等于每次运行都会调用所有函数。

| 分类 | 函数 | 触发条件或作用 |
|---|---|---|
| 默认配置 | displayDiffractionDefaults、coeLayerDefaults、coeLayerStackDefaults、coeSourceDefaults | 初始化配置；每次运行需要，不是一次性工具 |
| 文件读取与膜堆构建 | loadCoeProjectConfig、loadCoeSpectrum、resolveCoeImageFile、coeMaterialSpectrum、prepareCoeLayerStack | 按项目配置读取实际图样和材料；图层/光谱方式决定具体分支 |
| 照明与传播 | buildCoeIllumination、simulateCoeHeightStates、simulateDisplayDiffraction、propagateCoeReceiver、fourierIntegralCoe、reflectCoeFilmStack | 当前统一计算链；零高度只计算一个高度状态 |
| 指标与采样提示 | computeCoeAngularMetrics、computeCoeDiffractionSpots、reportCoeSampling | 角域、光斑、能量与采样评估 |
| 独立采样基准 | benchmarkCoeSampling、benchmarkCoePeriods | 主脚本runSamplingBenchmark=true时运行，mode="all"隔离五类因素，mode="periods"评估每轴Pitch数 |
| 绘图样式 | formatCumulativeAxes、formatFarFieldAxes、markDiffractionSpots、styleCoeFigure | 当前图形显示，标记开关只改变显示内容 |
| 条件色彩处理 | synthesizeCoeColor、composeCoeReflectionRgb | 前者用于D65/custom/monochromatic；后者仅RGB代理模式 |
| 历史批量接口兼容 | customPattern、normalIncidenceReflection | 主脚本不调用该旧几何/单界面分支；仅传cfg的历史runCOEReflection接口和兼容测试仍可能使用，当前不能算完全无用函数 |

以下5个文件已移至 `archive/input_generation/`，不再放在当前private目录：

| 归档文件 | 原用途 | 当前状态 |
|---|---|---|
| createCoeDiamondGeometry.m | 程序生成固定Diamond几何 | 已有用户图片，不在主流程调用 |
| createCoeImageExamples.m | 导出PNG/BMP/JPG示例输入 | 开发阶段数据准备工具，运行仿真不需要 |
| createCoeHeightExample.m | 生成可选二值高度条纹 | 演示工具；当前仅加载用户提供的heightMap |
| coeStackDefaults.m | 旧Anode/BM/CF固定结构默认值 | 已由coeLayerStackDefaults替代，归档前仅被旧工具和旧测试使用 |
| prepareCoeStack.m | 旧固定COE图样与膜堆入口 | 已由prepareCoeLayerStack替代，旧测试已迁移至统一图层接口 |

归档保存原始实现以便查阅，不加入MATLAB路径、不由主脚本/测试/打包入口调用。归档函数保留其原依赖，不承诺单独可运行；如需恢复，应在隔离副本中结合Git历史恢复相应上下文，不能直接将归档目录加入当前工程搜索路径。当前`src/`与`Input/`即可独立运行，不需要archive目录。

测试中的小尺寸图片和高度矩阵属于临时校验夹具，不是面向用户的输入生成程序；它们只在测试临时目录中创建，目的是验证实际文件读取、角域积分与回归行为，不覆盖用户输入。

### 程序责任与后续维护

| 文件或目录 | 用途 |
|---|---|
| src/COEReflection_main.m | 输入、配置、两状态计算、展示和保存 |
| src/coeReflectionDefaults.m、src/runCOEReflection.m | 保留兼容的公共配置/批量计算接口 |
| src/private/loadCoeProjectConfig.m | 项目文本配置解析 |
| src/private/prepareCoeLayerStack.m、reflectCoeFilmStack.m | 图样校验、区域查找表和相干膜堆 |
| src/private/simulateDisplayDiffraction.m | 复振幅、接收面积分与基础指标；独立FFT能量诊断 |
| src/private/fourierIntegralCoe.m、propagateCoeReceiver.m | 二次相位卷积傅里叶积分、传播系数、像素内强度积分 |
| src/private/loadCoeSpectrum.m、synthesizeCoeColor.m | D65、CIE积分和色度显示 |
| src/private/computeCoeAngularMetrics.m、computeCoeDiffractionSpots.m | 角域与光斑积分 |
| src/private/format*.m、markDiffractionSpots.m、styleCoeFigure.m | 图形样式和辅助标记 |
| Input/Colorimetry/ | 原始CIE数据、元数据、许可与校验值 |
| AGENTS.md | 本模块维护规则 |
| dist/ | 历史打包工具与旧快照，继续工作以当前源码为准 |
| archive/input_generation/ | 已退役的几何/图片/高度生成和固定结构入口，仅作历史源码存档 |

src 顶层仅保留主脚本与两个兼容公共接口。新增内部函数直接放在 private，使用小驼峰命名，文件名与首函数名一致；用户主流程和必要中文注释留在主脚本。现有数值算法回归由 tests/private/addCoeAlgorithmFixture.m 在临时副本运行，生产目录的私有解析由独立模块测试覆盖。

继续开发时优先接入实测图样与色散数据，核对CF剖面/平坦化和吸光底参数；可再增加手动ROI/阶次跟踪、角度截线与色度分布、有限光源/入射角/偏振/接收孔径。标准反射雾度需单独实现相应测量几何与校正。用户说明统一维护在本 README，不再按功能拆出多个平行说明文件。

### 后续提升计划：近场传播与人眼观察（暂缓）

当前已支持平面波、有限距离/尺寸点源及近轴TEM00光源，接收端仍为局部相干膜堆与标量近轴夫琅禾费角分布。当前阶段先完成光谱、空间及有限源采样收敛；近场传播和人眼成像暂缓，不增加相机模式。

后续需要定量对应实际观察光晕时，再分步扩展：

1. 有限距离点光源照明已实现球面入射波、光源位置/尺寸及非相干圆盘积分；后续结合实际5–50 cm光源距离核对照明条件和采样。
2. COE到观察位置的传播：增加Fresnel或角谱传播选项，并与现有夫琅禾费结果核对远场极限。
3. 人眼接收与成像：加入有限瞳孔、眼内透镜及对焦，关注眼睛距COE约20–50 cm时的反射衍射光晕。
4. 实测对照：统一光源光谱、观察几何和显示/曝光条件，比较光晕角宽、色彩及能量分布。

扩展前需重新确认实际光源、观察和对焦条件；当前归一化dB图用于比较结构衍射特性，不声明为实际视网膜亮度图。

批量接口还支持 `runCOEReflection(cfg,stack,inputDir)`，按与主脚本相同的列表构建膜堆；不自动加载config或执行全部后处理。只传cfg的历史接口仍保留。

```matlab
addpath('topic_coe_reflection/src');
cfg = coeReflectionDefaults();
cfg.complexCoefficient = 0.5; % 已知输入复振幅的演示，非实际材料参数
r = runCOEReflection(cfg);
disp(r.metrics);
```

<a id="faq"></a>

## 常见问题

- 默认只有一行图/一条曲线：heightMap全零时自动采用单状态模式，省去重复计算；非零高度才显示两状态对比。
- 零级降低但光晕增强：纯相位会重新分配能量，不一定降低总反射率；应结合观察角域、累计曲线和非零级能量判断。
- 找不到光斑：检查阈值、零级排除半径和积分半角，确认候选完整位于输出范围；没有候选时会输出只有表头的 CSV。
- CSV 中没有某些弱光斑：受到相对峰阈值、最小间距、最多光斑数及有限窗口筛选。这些选定光斑不代表所有非零级反射能量。
- 窗口能量不足 98%：中文警告逐项列出低于阈值的波长、窗口内/总反射比例、窗口外/总反射比例、窗口内/入射比例及当前观察角域。98%是提醒阈值，不是计算精度；窗口外能量不代表数值丢失。查看 `CapturedFraction`，按需扩大角度范围并相应提高 `outputSize` 保持角度分辨率；近轴模型限定半范围不超过10°。显示窗口外能量不会出现在光斑统计中。
- MATLAB 旧版本函数不支持：以已验证的 R2024a 环境为基准，不直接删除输入检查；先核对版本及语法兼容性。
- 中文字体异常：修改 `displayStyle.fontName` 为新电脑已安装的中文字体。
- 修改后结果异常：先关闭高度、保留默认参数运行验证；检查 m/nm、全角/半角、矩阵行列和波长顺序。

- 改主脚本参数却不生效：当前config会覆盖默认值，请改项目配置，或显式关闭useProjectConfig。
- 图片找不到或同名多格式：核对项目目录和文件名，重复格式时明确扩展名；所有层须同尺寸同原点。
- 同层交叠报错：核对layerLevels及invert，RGB CF应在同一层且掩膜互斥。不要为了绕过错误把三色配置成上下三层。
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
