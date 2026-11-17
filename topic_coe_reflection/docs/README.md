# COE 反射衍射：使用与交接说明

更新日期：2026-11-17。本文件是本模块的统一使用说明，包含项目配置、实际图样输入、物理模型、结果解释与开发交接。文中 `src/`、`Input/`、`Output/` 均相对于 `topic_coe_reflection/`。

正式仿真按图层列表读取用户提供的任意命名图样（默认阳极、BM、RGB CF），不调整像素形状。主入口是 `src/COEReflection_main.m`；选定项目后，自动读取该项目的 `config.txt`，无需反复修改 MATLAB 中的常用参数。额外高度默认 **0 nm**。

新增 `Input/F192Base/` 后，当前共25个项目，F192Base沿用下述统一运行参数；下文24项目/74图像的审计记录指新增前的快照。F192Base包含独立六层PNG、config与Materials，对应CAD为 `Input/Cad/F192/(20261103)F192Base.dwg`；同时更新COEB3的DWG，CAD备份与锁文件不纳入版本。

F192、COEB3及F192Base的Materials已保存新版BM/CF_R/CF_G/CF_B波长相关n/k CSV；F192/Materials同时保留四份原始TXT供追溯。数据格式为波长nm、n、k，CSV不含表头。三个项目当前 `nkFile=""`，运行时仍使用config中的n/k常数；归档CSV不会自动启用材料谱，启用路径和膜厚说明见下文材料章节。此处仅记录文件来源，不将未经核实的数据标注为实测。

**当前统一运行参数（2026-11-17，覆盖下文历史阶段参数）：** 全部25个Input config（含新增F192Base）已统一为4×4 Pitch、光栅→接收面25cm、出射角±5°、D65 380:10:780nm共41波长、平面波入射、无独立人眼入瞳。单波长内部空间相干，不同波长叠加强度；最终用CIE 1931 2° XYZ积分转sRGB。各项目原始图片、Pitch、层序、n/k、厚度及反射系数保留；当前验证继续运行case03_PlaneWave_Pitch4，通用主入口默认case01。下文中的旧周期数、RGB三波长和观察模型是历史对照记录，不代表现行config。

统一计算采样为64点/Pitch，角网格1801×1801。若原图片覆盖大于4×4 Pitch，config显式开启 `stack.geometry.cropToFit=true`，在降采样后的整Pitch边界上选中心区域；若原单元小于4×4，继续完整平铺。`inputSamplesPerPeriod`记录原生采样，保持物理Pitch；少数原生16点/Pitch的示例显式开启 `pixelReplication=true`，只进行整数倍像素块复制。两选项在代码默认均false，原图片文件不改写。`layerGeometry`记录inputPeriods、selectedPeriods、cropBoundsPixels、repeatCounts和pixelReplicationFactor。全24项目配置与图层几何已检查，74个输入图像哈希未改变；配置审计表在 `Output/ConfigAudit_380_780_10nm/config_audit.csv`。

全光谱不再把620/550/450nm直接指定为三个RGB通道；使用D65功率谱、10nm梯形积分权重（两端半权重）及CIE配色函数合成XYZ，再转sRGB。保留线性曝光sRGB与增强显示sRGB，增强仅按亮度Y统一缩放，参数toneStrength=60、luminanceGamma=0.72；旧legacy三通道显示参数只用于历史RGB模式。接收面毫米坐标由X=zθ得到；z=25cm、±5°对应X/Y范围约±21.8166mm、窗口全宽43.6332mm。仍沿用当前夫琅禾费近轴模型；380nm短波端的常用远场尺度约0.2937m，高于指定0.25m，程序如实报告此近似条件，不自动更改用户距离。

旧配置已冻结为 `tests/fixtures/coe_legacy_configs/*.legacy.txt`，旧分支的独立复制回归使用这些固定参数，避免随着当前用户preset更新而悄然改变测试含义；新增当前preset回归直接读取全部现行config。历史结果和历史参数不批量改写。

当前case03全谱结果已实际完成41个波长、1801×1801网格计算，最大Parseval误差约2.50×10⁻¹⁶。所有带坐标图的总标题通过 `plotLabelProject` 添加项目名称，原波长/周期/距离标题保留；无坐标RGB数据图保持干净像素。输出根目录的 `d65_rgb_receiver_mm_titled.png` 为带项目名及参数的接收面图，`d65_rgb_angle_titled.png` 为角度图，`d65_rgb_linear.png` 为线性曝光sRGB。配置/数值/历史分支独立运行共113项检查通过，加入项目标题后另通过当前全谱独立运行检查；记录在 `Output/ConfigAudit_380_780_10nm/final_regression.mat`。

可编辑的详细流程图：[topic_coe_reflection_flow.drawio](topic_coe_reflection_flow.drawio)。用 draw.io 打开后，依次查看底部四个页签：总览与项目数据流、图层到局部膜堆、照明与远场传播、指标色度与输出。图中标注实际函数名、主要公式、数组尺寸、能量归一化分母和当前默认 `Case1_Pitch64` 参数实例；内容依据 2026-11-16 工作区源码整理。

## 目录

- [快速开始与项目切换](#quick-start)
- [距离缩放观察视图](#distance-scaling)
- [有限光源与人眼入瞳](#finite-eye)
- [逐源点复场与薄透镜像面](#coherent-eye)
- [收敛性结论、当前选型与验证记录](#convergence-decisions)
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

主脚本默认 `showProgress=true`，显示中文仿真进度窗口：包含读图/膜堆计算、逐波长与源点传播、非零高度的两状态、指标统计、绘图及保存；显示当前波长、源点序号和已用时间。百分比按任务阶段分配，不表示精确耗时占比或剩余时间。完成或报错后自动关闭，手动关闭窗口只隐藏进度，计算继续。可在主脚本或项目 `config.txt` 设置 `showProgress=false` 关闭；无 MATLAB 桌面时自动跳过窗口。采样基准与批量公共接口沿用原有输出，不弹出该窗口。进度回调不写入 `cfg`、原始 XYZ 或物理结果。

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
  Case1_Pitch64/      Case1原图4×4平铺，4096×4096，完整64×64 Pitch
  Case1_Pitch128/     Case1原图8×8平铺，8192×8192，完整128×128 Pitch
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

Case1/2/3原始图片均为 **1024×1024**，Pitch=50 μm，输入单元为16×16 Pitch，每Pitch 64点、单元800×800 μm。当前config显式开启repeatToFit，按4×4平铺到64×64 Pitch，计算物面为4096×4096点、3.2×3.2 mm。参考图的实线 PDL 作为阳极，虚线是 Mask，不当作 BM。R/B 在网格上交错，G 位于四邻中心并交替倾斜；每个2×2晶胞有2R、4G、2B。BM 开口和 CF 按先前约定分别外扩5 μm和7 μm。参考图未标注尺寸，示例只匹配布局及轮廓特征，不代表产品尺寸或实测材料。

新增项目时复制一个 Case 文件夹并改名，替换六张图片及配置，再修改 `projectName` 即可（支持中文名称）。物理尺度、采样和材料由项目配置指定，不由文件夹名称推断；程序不缩放图片；只有显式开启stack.geometry.repeatToFit时才按下述规则平铺。`projectName=""` 使用旧的 Input/Output 根目录，启用配置时仍需在 Input 根目录提供 config.txt。

默认使用 D65 380–780 nm、5 nm 步长（81波长），输出1201×1201、角域±6°、接收距离10 m。有限圆盘光源默认每轴5点，光源类型仍由项目选择，默认plane。当前反射分支采用指定频率网格的二次相位卷积傅里叶积分；`paddingFactor` 保留兼容，不再控制反射输出采样。`cfg.heightMap=0`；只计算和展示一组反射结果，膜层自身的相位始终保留。两个单图层校验项目保留小尺寸输出，光谱同样使用5 nm。

<a id="distance-scaling"></a>

## 距离缩放观察视图

本支路用于验证有限距离点光源的**几何缩放假设**，不是完整人眼/相机仿真，也不是对实拍符合性的证明。光源距离为 `cfg.source.distance=a`，观察距离为 `cfg.observer.distance=b`，均为到屏幕的法向距离（m）。近轴、居中条件下，内禀衍射角q映射为观察角theta=s*q，s=a/(a+b)；屏幕投影位置为x=b*theta，故标称周期p的第m级角位置为s*m*lambda/p。

### 运行距离对比

在 `src/COEReflection_main.m` 选择 `projectName="case03_DistanceScaling"` 后运行。默认普通项目仍为case01。新项目只复制case03原始四圆Anode/BM图样，没有改变用户原项目；按已有显式平铺规则计算8×8 Pitch，保留250点/Pitch。配置为620/550/450nm三波长RGB代理、601×601输出、固定观察角域±1.5°，光源和观察距离基准均25cm。

```text
cfg.source.type = "point"
cfg.source.distance = 0.25
cfg.source.diameter = 0
cfg.source.integrationMethod = "direct"
cfg.observer.enabled = true
cfg.observer.distance = 0.25
cfg.observer.angleHalfRangeDeg = 1.5
cfg.observer.runBenchmark = true
cfg.observer.scanDistances = [0.15, 0.20, 0.25, 0.30, 0.35]
```

扫描先固定b、改变a，再固定a、改变b；输出位于 `Output/case03_DistanceScaling/DistanceScaling/`：

- `distance_scaling_comparison.png`：两排图分别移动光源和观察者，全部采用相同观察角坐标和一次共同曝光；不会分别放大图片抵消距离差异。
- `distance_scaling_peak_validation.png`、`distance_scaling_peaks.csv`：550nm标称Pitch第一阶候选峰的理论位置与能量图读数。候选区为预测位置的0.7至1.3倍；仅适用于输入确有该级次的情况，不根据文件名或图像猜测最小晶胞，不代表任意图样都有该峰。
- `distance_scaling_results.mat`、`config_snapshot.txt`：各距离的原始逐波长参考能量、角坐标、缩放系数、图样/材料及配置记录。

`runBenchmark=false` 时执行普通单项目计算，并在原有远场结果之外增加独立观察视图；输出为 `*_distance_view.png`，原始数组保存于结果MAT的 `result.distanceView`。批量API同样支持：

```matlab
cfg = coeReflectionDefaults();
cfg.source.type = "point";
cfg.source.distance = 0.25;
cfg.observer.enabled = true;
cfg.observer.distance = 0.25;
result = runCOEReflection(cfg); % 也可传入原有stack和inputDir
view = result.distanceView;
```

### 数值与物理边界

`coreDistanceView` 用相同图样、复反射系数和明确的plane照明足迹，重新计算平面波参考在q=theta/s处的傅里叶积分，再将坐标映射到观察角。它在每个波长的原始能量上运算，不缩放已增强的RGB图，不重复叠加球面入射相位。`cfg.source.wavefront` 用于原有远场支路，不改变这个几何参考；照明足迹由source.apertureShape/rectangleSize/circleDiameter定义，须在距离比较中保持固定，且不模拟点源照明衰减。

对应角像素的能量不变；角密度满足J_view(theta)=J_ref(theta/s)/s²，`view.spectralAngularDensity` 已使用映射后的角像素面积。固定观察窗截断能量时不重新归一化。`capturedReferenceFraction`和`referenceThroughput`都是平面波参考量；本支路不输出眼睛接收功率或W/m²。`screenCoordinate=b*coordinate`是视线投影到屏幕平面的坐标，不是远场接收面光斑位置。

`cfg.propagationDistance`仍是原有夫琅禾费接收面的数值距离，与眼睛距离分开，不能把它改成25cm来替代此模型。原有远场能量、雾度代理指标及辐照度不会被几何视图覆盖。几何支路要求居中point及direct，已支持下节的有限圆盘源和独立人眼入瞳近似；偏轴和激光仍拒绝。包含卷积补边后的参考半角必须不超过10°并通过物面Nyquist检查。尚未包含完整相干眼内成像、像差、离焦、光源封装光学和实际屏幕边缘对视野的裁切。

解析相位斜坡测试通过公共API定位实际峰值，验证不同a/b的角位置、屏幕投影位置、密度Jacobian、窗口截断不补能量、原远场结果不变及独立复制解析。测试以数值积分与独立解析峰位比较；距离模型符合实拍的程度仍需测量。

2026-11-16原生四圆示例实算：550nm、p=59.06μm、观察角网格0.005°。固定b=25cm，a=15→35cm时，候选峰读数0.200→0.310°；固定a=25cm，b=15→35cm时，读数0.330→0.220°。十组理论/读数最大绝对差0.003482°，小于一个观察角网格间距；有限阵列包络也可能使局部峰偏离无限周期级次。±1.5°固定窗口捕获约91.66%至95.22%的550nm参考反射能量，未将截断能量补回。对比图及原始CSV位于上述输出目录。

本轮全仓库104项MATLAB测试全部通过，其中新增12项公共API数值/输入/复制检查，以及3项独立复制模块的距离扫描、RGB主脚本和光谱主脚本检查。完整日志保存在上述输出目录的 `regression_tests.log` 和 `regression_tests.mat`。

<a id="finite-eye"></a>

## 有限光源与人眼入瞳

在主脚本选择 `projectName="case03_FiniteEye"`，运行同一入口即可比较光源距离和观察距离各15、20、25、30、35cm的结果。光源发光直径固定5mm，人眼入瞳直径固定3mm。以下参数在该项目config中独立保存：

```text
cfg.source.diameter = 0.005
cfg.source.distance = 0.25
cfg.observer.pupilDiameter = 0.003
cfg.observer.distance = 0.25
cfg.observer.pupilModel = "airyIntensity"
cfg.observer.airyRadiusFactor = 20
cfg.observer.referenceMethod = "fftBins"
cfg.observer.referencePaddingFactor = 1
```

`observer.pupilDiameter`是眼睛的有效入瞳全直径，不修改屏幕计算面上的 `cfg.pupilRadius`。`pupilModel="geometryOnly"`仅统计几何贡献区域；`"airyIntensity"`额外对逐波长角图使用理想圆形入瞳的非相干强度PSF。默认geometryOnly保留此前点源角图；人眼入瞳默认3mm。直径0表示此前点观察者的数学极限，并非物理上零孔径仍可接收光。

### 图样计算与适用范围

先计算theta=a/(a+b)*q的点源参考，再按均匀非相干圆盘源进行逐波长强度平均。源的观察角全直径为Ds/(a+b)，不是Ds/a；源角核积分为1。独立源点的强度叠加与相干场求和需要区分，参见[Ansys非相干线源示例](https://optics.ansys.com/hc/en-us/articles/42661828176147-Simulation-of-Young-s-interference-experiment-via-geometric-ray-tracing-in-OpticStudio)。

本观察视图的源积分采用角网格上的圆盘概率核，`source.samplesPerAxis`不控制这一步；它仍用于原有direct远场支路的离散源点求和。此视图应通过观察角网格、参考FFT分箱和角核采样检查精度。

理想圆孔的强度PSF正比于[2J1(v)/v]²，v=pi*Dp*theta/lambda，第一零点半角约1.21967*lambda/Dp。3mm入瞳在550nm时约0.0128°。此处**显式采用移不变、非相干角图强度PSF近似**；它不是一般部分相干光的完整眼内传播，不把瞳孔内各点视作独立发光源，也未模拟调焦和像差。一般相干成像需要保留复振幅传播，参见[Ansys相干/部分相干成像说明](https://ansyshelp.ansys.com/public/views/secured/zemax/v26103/en/opticstudio_user_guide/OpticStudio_Help/topics/Partially_Coherent_Image_Analysis.html)。该近似是否符合实拍仍需测量。

卷积先按源核与入瞳核总支撑半径向外补边，在逐波长参考能量上做FFT线性卷积，最后裁切固定观察角窗口；不对RGB增强图做模糊，不补回窗口外能量。圆盘源核使用子像素面积采样；Airy核按第一零点的20倍截断，离散积分校准到解析包围能量1-J0(v_edge)²-J1(v_edge)²，约0.992，而不是将截尾核归一为1。核外约0.8%功率保留为明确的截尾损失。`view.kernels`保存角核、第一零点、核积分及补边尺寸，`convolutionMassError`检查完整线性卷积的功率记账。此分支输出仍是参考能量分布，不是眼睛实际接收瓦数或视网膜辐照度。

`observer.referenceMethod="quadrature"`保留原来的指定频率傅里叶积分，仍是旧点源缩放项目的默认算法。新大计算面项目显式选择`"fftBins"`：按FFT原生角像素的面积重叠分配能量，再执行宽面源卷积；原生角像素内部采用常密度近似。`referencePaddingFactor`为FFT零填充整数倍，可用1/2做数值敏感性检查。该近似适合当前宽面源定性图样，不用于替代原有精细点源峰形验证；未经检查不能把零填充当作真实图样采样细化。两种算法保留独立记录，原有远场支路默认不变。

### 几何贡献区域

居中、近轴平面镜模型下，源点u、入瞳点v对应的屏幕反射点为x=(b*u+a*v)/(a+b)。对两个圆盘进行投影后，光线对位置分布是两个圆盘的卷积；支撑区域直径为：

```text
L = (b*Ds + a*Dp)/(a+b)
N_nominal = L / cfg.pitch
```

`view.participation`独立保存源投影直径、入瞳投影直径、支撑直径、面积、标称Pitch数以及401×401诊断网格上的光线对密度。`rayWeightInModel`为几何光线对落在当前计算面与照明足迹交集内的比例；`modeledArea`为该交集的诊断网格面积近似。**该区域只描述镜面零级，不是所有衍射级次的统一硬孔径；标称Pitch数不是相干周期数。** 程序不据此裁切反射场，也不将它当作衍射能量归一化分母。

两距离均25cm时L=4mm，面积约12.57mm²，跨越约67.73个59.06μm标称Pitch。固定观察距离25cm，光源距离15→35cm时L从4.25缩至3.833mm；固定光源25cm，观察距离15→35cm时L从3.75增至4.167mm。图样的角尺度与这个屏幕区域宽度是不同量。

### 新示例和输出

`case03_FiniteEye`原样复制case03的两张500×500图片，不改写旧项目。为覆盖最大约4.25mm的几何区域，config显式平铺到80×80 Pitch，即4.7248×4.7248mm；显式用inputSamplesPerPeriod=250、samplesPerPeriod=32做共同像素中心最近邻降采样，计算网格2560×2560。参考照明apertureShape=none表示该完整计算面参与计算。此采样用于当前距离和有限源平均的定性比较，不宣称原生轮廓精度收敛；源核与眼PSF也不反过来改变用户图片或自动扩展任何其它项目。

在 `Output/case03_FiniteEye/DistanceScaling/` 查看：

- `distance_scaling_comparison.png`：距离扫描的最终RGB角图，全图共同曝光和固定坐标。
- `observer_source_pupil_stages.png`：25cm基准下点源、5mm源、3mm入瞳近似三阶段的共同曝光RGB和550nm截线。
- `observer_participation_regions.png`：镜面几何区域及实际输入边界，两排分别扫描光源和观察距离。
- `observer_geometry_summary.csv`：区域直径、面积、标称Pitch数、计算面覆盖比例、源角宽度和候选峰读数；MAT保留全部逐波长数组及核。

有限源角宽度达到预测级次间距的0.8倍时，`PeakComparable=false`，最终图的 `MeasuredPeakDeg`/`PeakErrorDeg`写NaN；这个保守可比较性判据不是通用分辨率标准。点源参考峰另存 `PointReferencePeakDeg`。不能从源平均后的混合亮区强行指认原来的单级峰；源导致点阵被洗平也不能用减小源直径或逐图曝光来掩盖。

2026-11-17实算：10组均为几何区域完全落入4.7248mm计算面的情况，rayWeightInModel=1。5mm均匀源明显将点源点阵展宽、叠合；基准观察角全直径0.57296°，而550nm第一阶间距约0.26679°，因此最终图的单级峰均标记为不可直接比较。3mm入瞳的第一零点半角0.01281°，在当前近似下进一步平滑的影响小于源尺寸平均。该结果尚未匹配实拍清晰点阵，也不能由此唯一判定实拍的光源相干性或封装光学。

同一32点/Pitch几何、相同角核和固定观察窗口下，FFT原生分箱与原指定频率积分的10组×3波长最终能量图，最大相对L1差异1.08145%，最大窗口积分能量相对差0.002287%。对照数据为项目Output中的 `quadrature_reference.mat`、`fft_vs_quadrature.mat` 和 `fft_vs_quadrature_relative_L1.csv`。这只验证分箱近似的影响，不是对32点/Pitch轮廓精度或实拍模型的验证。

本轮全仓库117项MATLAB测试全部通过；包含公共接口下的独立卷积对照、源角核二阶矩、观察窗外峰卷入、源/入瞳区域公式及光线对密度积分、Airy解析包围能量、零尺寸退化、入瞳与屏幕孔径独立、FFT分箱近似、复制模块运行，以及原有COE/UDC等分支回归。日志为本项目DistanceScaling输出目录的 `regression_tests.log`、`regression_tests.mat`。

### 更宽的观察角域

选择 `projectName="case03_FiniteEye_WideAngle"` 可运行光源/观察距离均25cm的宽角基准，仍为5mm均匀源及3mm入瞳。每轴观察半角扩大到4°，总宽8°；1601×1601输出保留0.005°角步长。物面采样显式提高为64点/Pitch（5120×5120），以容纳卷积补边后的±9.16°参考角域，并保留原80×80 Pitch和4.7248mm物理尺寸。图样图片未改写，旧±1.5°项目和结果保留。

`Output/case03_FiniteEye_WideAngle/wide_angle_fov_comparison.png`对同一次±4°计算裁切显示±1.5°、±3°、±4°，使用相同曝光和角度步长；它展示更多外围角域，不把中心图样拉伸。`wide_angle_4deg_rgb.png`为完整宽角图，`wide_angle_summary.csv`记录各窗口的参考捕获能量。绘图脚本在该输出目录的 `render_wide_angle_views.m`，数值结果和配置快照在 `DistanceScaling/`。直接运行主脚本会生成标准宽角结果图；三视场对照由该绘图脚本读取同一MAT结果生成。

本宽角项目只扫描[0.25]基准距离。若恢复15至35cm全距离扫描，最小a/(a+b)下的补边参考角域会更大；当前参数建议观察半角最多3°，继续遵守10°参考角与物面Nyquist约束，不能静默解除限制。宽角结果仍沿用非相干角图/理想入瞳近似，未验证实际屏幕边界接收或完整人眼成像。

### 小光源、减少周期与高光压缩

选择 `projectName="case03_FiniteEye_SmallSource"`：光源直径从5mm改为0.5mm，计算阵列从80×80减为16×16 Pitch；入瞳保持3mm，光源和观察距离均25cm。64点/Pitch、1024×1024计算网格对应0.94496mm物面。观察角域扩为每轴±4.5°，1801×1801输出，仍为0.005°角步长。显式FFT零填充倍数提高到4以检查较少周期下的峰形。输入图片直接复制，单元形状不调整。

这里减少的是实际建模的阵列范围，不是将相干周期数指定为16。镜面几何支撑直径为1.75mm（约29.63个标称Pitch），大于0.94496mm计算面；程序如实报告交集面积和落入计算面的光线对比例，不自动补齐被用户刻意缩减的区域。

RGB代理新增可选显示模式，旧项目默认仍为legacy：

```text
rgbOptions.mode = "luminance"
rgbOptions.clipQuantile = 1
rgbOptions.alpha = 10000
rgbOptions.gamma = 0.8
rgbOptions.peakLevel = 0.7
```

新模式用0.2126R+0.7152G+0.0722B作为RGB代理的显示亮度，以同一参考值做单调对数映射，每个像素的三个通道共用缩放，再对RGB向量整体压缩以满足最大通道0.7。它提高弱信号相对于高光的可见度，保留代理通道比例，不按中心坐标设置暗斑或遮罩；不宣称三波长代理具有严格CIE色度。原始spectralEnergy、光斑位置、积分能量均不修改，连续光谱的XYZ流程不受此RGB显示选项影响。

标准增强图保存为 `Output/case03_FiniteEye_SmallSource/DistanceScaling/observer_baseline_rgb.png`；同一目录保存原始能量和显示参数。输出根目录的 `small_source_contrast_comparison.png`、`small_source_legacy.png`、`small_source_enhanced.png`由同一份新物理结果比较两种显示曲线，不把显示差异解释为物理能量变化。绘图脚本为同目录的 `render_contrast_comparison.m`。

本项目采样检查：相对于4倍FFT零填充，不填充的全角图相对L1差异约20%至23%，因此最终配置使用4倍。中心±1°以相同0.005°网格与指定频率直接积分比较，三波长最大相对L1差异小于0.6%；此检查不等于对全角域、连续轮廓或实拍模型的精度保证。记录为 `sampling_check.mat`。显示增强与旧项目RGB/光谱流程共11项相关检查通过，包括独立复制主脚本运行、0.7通道上限、RGB向量比例在PNG量化误差内保持；记录为 `display_tests.log`、`display_tests.mat`。本轮未重复无关物理分支的全仓库测试。

### 2mm光源与48×48周期

选择 `projectName="case03_FiniteEye_Source2mm_Pitch48"`：光源直径2mm，计算阵列48×48 Pitch，观察角域每轴±4.5°；入瞳3mm、光源和观察距离均25cm。64点/Pitch，3072×3072计算网格，物面2.83488×2.83488mm。使用2倍FFT零填充，当前显示采用早期case03的legacy模式。卷积补边后的参考半角约9.81°，仍在10°限制内。

所有正式观察角图、阶段图和峰位图的标题自动读取cfg.periods，标注实际计算周期数，例如“48×48 Pitch”。几何区域图另外标注其支撑直径折合的Pitch数：本组为2.5mm、约42.33个标称Pitch；该值与实际输入的48×48计算阵列区分。

标准结果在 `Output/case03_FiniteEye_Source2mm_Pitch48/DistanceScaling/`，带参数和周期标题的单张增强图为输出根目录的 `source2mm_pitch48_titled.png`，无坐标原图为 `source2mm_pitch48_rgb.png`；绘图脚本为 `render_titled_view.m`，参数表为 `parameters.csv`。

本组最初执行检查确认光源2mm、周期[48,48]、半角4.5°及标题文本一致，旧图使用0.7高光上限。中心±0.75°与同角网格的指定频率直接积分比较，三波长最大相对L1差异约0.25%；局部检查记录于 `sampling_check.mat`，不代表全角域、轮廓采样或实拍模型均已验证。

对照实拍后，本项目曾比较温和显示 `alpha=300, gamma=1, peakLevel=1, mode="luminance", clipQuantile=1`。旧设置10000/0.8/0.7同时大幅提升弱信号、限制中心高光，使中心与外围的亮度差过小；温和配置减弱对数增强、取消额外伽马提亮并恢复高光上限。该选择仅用于定性比较亮暗层次，未经相机曝光或响应标定；仍不能靠显示曲线校正光斑尺寸、位置或三波长代理颜色。

温和显示直接复用同一份原始能量结果，不重算传播、不改写历史MAT及其mapping记录。对应图为输出根目录的 `source2mm_pitch48_mild_titled.png` 与 `source2mm_pitch48_mild_rgb.png`，脚本为 `render_mild_view.m`；独立的 `mild_display.mat` 记录显示参数、引用数据路径和检查结果。`milder_mapping_comparison.png` 比较alpha=10000/1000/300/100，上排完整±4.5°，下排裁看±1.5°，裁看不提高计算分辨率。历史标准图及 `render_titled_view.m` 保留原强压缩版本。2°至4°环带的显示平均亮度由约0.091降至0.00784，中心最大通道由0.7恢复到1；这些是显示量，不能解释为能量变化。

按用户选择，当前项目最终切换为早期case03的 `mode="legacy", clipQuantile=0.9985, alpha=60, gamma=0.72, peakLevel=1`。它对三个通道共用分位数参考值，但分别做截断、对数及伽马映射，因此可出现白色高光，也会改变通道比例，不能用于判断严格色度。光源2mm、48×48 Pitch、入瞳3mm、a=b=25cm和±4.5°等物理设置保持不变。当前重绘文件为 `source2mm_pitch48_legacy_titled.png`、`source2mm_pitch48_legacy_rgb.png`，脚本为 `render_legacy_view.m`，独立显示记录为 `legacy_display.mat`；原始传播MAT及此前显示版本保留。未来主入口运行本项目自动采用config中的legacy设置；当前视图的分位数参考由本视图信号确定，早期四图样对比使用全部图样共同参考，两者不保证相同曝光。

<a id="project-config"></a>

### 4×4 Pitch有限阵列与强度卷积对照

选择 `projectName="case03_FiniteEye_Source2mm_Pitch4"`：在2mm光源项目基础上，将 `cfg.periods` 改为[4,4]。Pitch59.06μm，物面边长0.23624mm，每Pitch64点，256×256网格；观察±4.5°，1801×1801输出。光源直径2mm、入瞳3mm、a=b=25cm、legacy的0.9985/60/0.72显示参数保持一致。FFT参考零填充16倍以细化小阵列宽峰的能量分箱；补边规则、能量归一约定沿用距离模型。48×48项目保留，主脚本默认仍为case01。

基本物理计算顺序是：对每个独立源点计算相干复场传播，再将各源点产生的强度按权重相加。只有源点变化可近似为同一强度图的角度平移、照明足迹和反射响应不随源点改变时，才可把该积分简化成角域强度卷积。本距离视图继续采用这一近轴平移近似；核直径为Ds/(a+b)，本组约0.22918°。这不是对复振幅卷积，也不意味着物面各点互不相干。入瞳Airy强度卷积是另加的非相干角图近似，并非完整相干眼成像。

扩展光源的源点强度平均本身已经包含其造成的条纹可见度降低，不能再把由同一光源推导的相干长度作为独立模糊核重复施加。本次4×4只作为用户指定的有限计算孔径敏感性对照，没有额外施加相干包络，也没有实现“大面积屏幕由4×4互不相干小块组成”的统计模型。若将4×4解释成真实有效相干范围，应另建立互强度或相干模态模型并校核，不能由本图直接断言。

镜面零级几何贡献区域仍为2.5mm跨径，明显大于0.23624mm计算面；相交面积和光线对权重由几何诊断独立报告。有限阵列图仅比较角分布形状，不代表整块屏幕亮度或接收功率。输出在 `Output/case03_FiniteEye_Source2mm_Pitch4/DistanceScaling/`，逐波长保留点源参考、光源平均和入瞳处理后的数组，可判断各阶段对渐变边缘的影响。

带标题的完整观察图为输出根目录 `source2mm_pitch4_titled.png`，三阶段形状图为 `pitch4_source_pupil_shapes.png`，由 `render_pitch4_views.m` 绘制。三阶段图各自使用legacy分位数参考，仅比较形状，不比较绝对亮度。光源平均后仍可见比48周期版本更平缓的亮斑边缘，入瞳阶段变化较小。中心±0.75°用同角度分箱的直接积分交叉检查：点源参考三波长最大相对L1差约1.25%，带补边执行光源与入瞳卷积后的最大差约0.059%；不代表全角域或真实模型精度。记录和复现脚本为 `sampling_check.mat`、`check_pitch4_sampling.m`。主计算卷积质量误差小于1e-15，几何光线对落入当前计算面的比例约0.0328，观察窗口捕获比例约0.9727（550nm）。

<a id="coherent-eye"></a>

### 逐源点复场、圆瞳与薄透镜像面

选择 `projectName="case03_CoherentEye_Pitch4"`，并使用 `cfg.observer.model="coherentEye"`。该分支从图片膜堆构建复反射系数，对每个非相干源点逐波长执行：有限距离入射复场 → 反射复场 → 距离b的Fresnel复场传播 → 瞳面乘圆孔及薄透镜相位 → 像面复场传播 → 模平方 → 按源权重求和。不同波长在取强度后合成三通道RGB代理。不会额外施加旧a/(a+b)缩放、光源强度卷积或Airy强度卷积，旧项目仍使用geometric默认分支。

当前输入保持2mm均匀圆盘源、3mm入瞳、a=b=25cm、4×4 Pitch、64点/Pitch和legacy显示。用户选择对焦镜中光源，`focusMode="virtualSource"`令对焦物距为a+b=50cm；另支持 `screen`（b）和 `distance`（显式focusDistance）。理想透镜到像面距离由 `retinaDistance=0.017` m定义，焦距满足1/f=1/v+1/L；此像距是模型坐标标定，不是对眼睛处方的测量。输出横纵坐标为像面位置/v的近轴角，不直接把衍射出射角当成人眼看到的像位置。

Fresnel复傅里叶积分保留dx*dy、1/(iλb)和相位；透镜相位与像面传播相位可解析合并为瞳面上的exp[-iπ|p|²/(λL)]，像面角强度密度等于该瞳面复场傅里叶积分的模平方除λ²。按像素中心及角面积Δθ²求积，不将输出窗口归一化为1；记录原反射功率、进入瞳孔的功率与像面窗口捕获功率。所有功率比例相对于有限受照计算面的入射参考功率，跨项目窗口大小不同时不能据此比较整屏绝对亮度。

为满足最短450nm波长的采样，完整±4.5°像面使用2401×2401点（0.00375°步长），圆瞳每轴1537点。程序分别检查瞳面复场/输出角域Nyquist和像面强度带宽；不足时明确报错，不静默产生周期复制图。输出的局部裁看使用相同图像和曝光，不提高真实分辨率。源采用 `source.quadrature="polarGauss"`：t=r²/R²上的7点Gauss与32个方位角，合计224个互不相干源点；面积权重总和1，不添加随机相位或额外相干包络。旧项目默认cartesian取样不变。

输出在 `Output/case03_CoherentEye_Pitch4/CoherentEye/`：`coherent_eye_results.mat`保留逐波长像面能量、源坐标及权重、瞳面辐照度、功率比例、对焦/透镜参数及显示映射；`coherent_eye_comparison.png`同时显示±4.5°与±0.25°裁看；`coherent_eye_rgb.png`为无坐标原图。额外高度全零时只算一个状态；非零时在单源复反射场上添加高度相位并分别计算两状态，不提前取模平方。

完整视野绝大部分为空时，全图0.9985分位数会把中心过曝。此项目用 `displayReferenceHalfAngleDeg=0.25`指定共同曝光参考ROI；只用该区域确定一个全图共用的参考值，不设置空间亮度遮罩，不更改能量。legacy的60/0.72与0.9985参数不变。`eyeDisplay.mapping`记录参考值及ROI，0表示完整图参考。首次全图曝光的图像另保留为 `coherent_eye_full_frame_exposure.png`，供区分曝光差异与物理差异。

本次完整计算每波长224个源点、零额外高度，用时约405秒。620/550/450nm进入瞳孔的功率/有限物面参考入射功率分别约0.03885/0.04291/0.05628；像面±4.5°窗口/瞳内功率分别约0.99989/0.99983/0.99984。4×4有限面下主要呈现中心小角域亮区，没有复现旧角图的大范围彩色点阵。这是应继续检验的几何范围假设，不能靠显示曲线修补。

这仍是法向局部反射、标量近轴、无像差薄透镜模型，不是矢量电磁仿真或完整人眼生理模型。尤其4×4 Pitch物面只有0.23624mm宽，会限制能够进入眼睛的衍射路径；不能把它理解成保留大屏幕几何范围、只把相干长度设为4个周期。若要模拟大屏幕上的部分相干，应保留相应物面范围并通过同一源积分或互强度模型处理，不能根据旧图外观再施加独立模糊。

源积分检查在450nm、4×4原图片、±0.3°和257点瞳面下对比7×32与9×48的圆盘积分，差异小于0.01%；记录在输出根目录 `validation.mat`。在相同像面角网格上，257/513点瞳面相对1537点参照的蓝光±0.3°图样相对L1差分别约0.19%/0.058%；复现脚本为 `check_pupil_sampling.m`。新增12项测试检验复傅里叶积分的相位、Cartesian及非等权Gauss源强度和、相位/反射率缩放、对焦位置与像反转符号、额外高度、圆盘积分矩、采样拒绝、固定曝光参考、黑反射面与旧API兼容；另1项独立复制模块运行main验证private解析和项目配置读取。既有距离/膜堆/光谱/基础衍射93项回归均通过，记录在 `regression_tests.mat`。新增核心文件Code Analyzer无提示。采样检查范围及精度不能代替真实光路校准。

### 退回纯平面波、取消人眼入瞳

当前对照项目 `case03_PlaneWave_Pitch4` 按用户要求采用 `source.type="plane", diameter=0, wavefront="plane"`，并设置 `observer.enabled=false, observer.pupilDiameter=0`。此处取消的是独立人眼入瞳和观察成像分支；屏幕计算面沿用默认 `cfg.pupilRadius=Inf`，仅保留用户指定的4×4 Pitch有限反射范围。config不将Inf写成非JSON字面量，也不用极大的有限数代替，因为有限光瞳须完整落在物面范围内。

用户已确认后续以本项目作为计算基准：平面波、4×4 Pitch、±10°、无独立入瞳及人眼成像、legacy显示；暂不继续圆盘源显示。日后讨论距离或相干性时，先与这个固定基准作对照，不自动恢复旧圆盘/几何缩放/入瞳分支。主脚本的通用默认case01仍保留；运行本基准时显式选择本项目。

当前图是远场角分布，不计算实际观察者位置：source.distance对plane不生效，observer.distance在observer.enabled=false时不生效；已有a/(a+b)模型也已关闭。cfg.propagationDistance=100m是本项目的远场接收面换算距离，不是25cm人眼观察距离。相同角度网格下，改变此接收距离不会改变归一化角分布；接收面线性坐标X≈zθ及光斑长度随z增大，像面每单位面积辐照度按1/z²变化。这些结论属于本模型的理想远场，不等于实际屏幕观察的距离效应已经被验证。

按用户要求，当前基准另开启 `cfg.receiverView.enabled=true`，用 `cfg.receiverView.distances=[0.25,0.5,1]` m输出25/50/100cm接收面的毫米尺寸。该列表中的z明确表示光栅→接收面，复用同一远场角分布，不重新传播、不更改光源→光栅距离a；原100m计算结果只是辐照度换算的参照。毫米坐标X=zθ；每个角像素在接收面的面积为(zΔθ)²，辐照度由参考结果乘(z_ref/z)²，角窗口内接收功率不变。本对照是远场假设下的坐标换算，不是对25cm处Fresnel近场分布的验证。

新增输出 `coe_reflection_receiver_rgb_mm.png` 显示列表中首个距离（默认25cm）的毫米坐标RGB图；`coe_reflection_receiver_distance_scaling.png` 在统一毫米坐标上下限上比较多个距离，避免各面板自动缩放掩盖尺寸变化。RGB图采用相同归一显示以比较尺寸，实际1/z²亮度变化另记录在 `coe_reflection_receiver_distance_scaling.csv`；MAT的 `displayData.receiverDistanceView` 保留毫米坐标、像素面积、辐照度缩放及汇总表。550nm标称一级位置x≈zλ/p在25/50/100cm分别为2.328/4.656/9.313mm；有限阵列的真实峰位可能受单元包络影响，此列是标称值。

每个620/550/450nm波长使用一个跨整个计算面的相干平面波，不同波长按强度求和。时间相干性降低应通过实际光谱及波长积分处理；给全场添加共同随机相位不会改变其衍射强度。空间相干性降低需要多个入射模式或互强度模型：若不同周期互不相干，跨周期干涉项消失，不能沿用当前完整阵列的场求和，也不能只把RGB显示做模糊。严格单方向、单频率理想平面波是空间相干的，需将其与宽带准直光或互不相干平面波模式的集合区分。

用户已明确“非相干平面波”指时间不相干、空间相干。后续按此含义处理：每个波长保留全物面的空间相位关系，不同波长叠加强度；不额外给各周期添加独立随机相位。真实时间相干程度由光谱宽度决定；当前三离散波长仅是颜色/光谱代理，不等价于已完成连续白光积分。严格平面波入射条件不变且接收面固定时，仅改变光源距离a不能产生X=zθ缩放；若要研究a的影响，需要有限距离照明波前或相应反射几何，单点球面/抛物波即可研究，不必引入圆盘源。

Pitch59.06μm，4×4阵列边长0.23624mm；每Pitch64点，256×256物面。按用户最新要求，输出由±4.5°扩为±10°，保留1801×1801角网格，角步长约0.01111°；450nm的λ/D尺度约0.10914°，仍有约9.8个显示间隔，像素内继续用3×3积分。物面最短波长Nyquist半角约13.97°，覆盖当前像素边界。此模型仍使用近轴角坐标，10°是当前实现的范围上限。入射为每个波长一份相干平面波，620/550/450nm以原D65相对权重按强度合成RGB代理。保留legacy的clipQuantile=0.9985、alpha=60、gamma=0.72。无圆盘源平均、源强度卷积、Airy强度卷积、瞳面/透镜成像或几何距离缩放。输出横纵坐标为光栅远场出射角，不能与上一项目的像面位置/像距坐标直接等同。source.distance对plane不生效，propagationDistance仅用于远场接收面坐标与辐照度换算。

标准结果位于 `Output/case03_PlaneWave_Pitch4/`；带参数标题的单张图为 `plane_wave_pitch4_titled.png`。旧圆盘源和人眼成像项目保留以供后续定位差异，主脚本默认case01保持不变。

### case01~case04 四图样平面波对比

用户指定四种理想图样，按截图从左到右命名：`case01`单圆、`case02`等面积椭圆（长短半轴比1.15）、`case03`四圆阵列、`case04`四个不同取向的等面积椭圆。各项目在 `Input/case01` 至 `Input/case04`，分别包含anode.png、BM.png和config.txt；图样是理想对照，不替代F192真实版图。由归档工具 `archive/input_generation/createFourPatternCases.m` 一次性建立，该工具拒绝覆盖既有项目；正式计算仅从图片读取，不提供形状调整参数。

截图中的毫米尺寸仅用于确定相对比例。本次统一缩放为中心间距Pitch=59.06 μm、圆半径Pitch/4=14.765 μm；椭圆半轴为r√1.15、r/√1.15，连续面积相同，二值采样后覆盖率会有微小差异。每张500×500图片覆盖2×2 Pitch（118.12 μm），保留250点/Pitch。case04的左上/右上/左下/右下旋转角采用60°/30°/−30°/−60°，是明确的示意假设，截图未提供精确角度。Anode图中白区反射率1、相位180°；BM与Anode同图，白色表示开口、invert=true，其余为零反射黑区。

主脚本当前选择 `projectName="case01"`、`runFourPatternComparison=false`，默认执行普通单项目的完整指标和绘图。手动设`runFourPatternComparison=true`时，统一读取并比较全部四组，然后返回。私有 `benchFourPatterns` 复用正式 `runCOEReflection` 传播。每组分别计算完整2×2 Pitch局部图样与8×8 Pitch阵列（完整图片每轴重复4次），均使用平面波、620/550/450 nm、±6°角域、1201×1201输出及相同D65三点权重。四组共八张RGB图共享同一个曝光参考与显示参数；灰度行仅比较局部图样，共享550 nm峰值参考、显示−60至0 dB。单圆和四圆的总反射面积不同，因此共同曝光图保留功率差异，没有按各自峰值归一化。

输出：各项目 `Output/case01` 至 `Output/case04` 保存 `local_rgb.png`、`array_rgb.png` 和 `pattern_metrics.csv`。统一 `Output/FourPatterns/` 保存 `four_patterns_comparison.png`（四列依次case01~04，四行依次输入图样、局部550 nm、局部RGB、重复阵列RGB）、`case03_case04_local_rgb.png`（局部四圆/四椭圆大图）、`four_patterns_summary.csv` 及 `four_patterns_results.mat`。MAT中的report保留各case局部/阵列的原始spectralEnergy、配置、入射面积与功率、重复数、局部有效反射图和共同映射；可用spectralEnergy乘cfg.spectralWeights重建RGB输入。图像仍是每单位总入射功率的角像素能量，未模拟光源大小、连续光谱或相机，不因视觉相似而宣称复现实拍。

2026-11-16实算：case01~04的总反射率分别为4.9124%、4.9100%、19.6496%、19.6304%，局部与重复阵列相同；最大Parseval绝对误差2.78×10⁻¹⁶，±6°窗口捕获约96.3%~97.4%的反射能量。局部单圆/单椭圆主要形成彩色环，四圆形成规则彩色光点，四个不同取向椭圆使环状结构破碎且光点强度随方向变化。重复后均有更锐利的点阵；case01/02此时也会产生阵列干涉，不能把“圆只产生环”推广到重复阵列。与实拍只作定性比较，不能由相似性唯一反推结构。新增四图样独立复制运行测试及原真实图样平面波复制测试均通过；测试用低分辨率临时夹具验证流程和共同曝光，不作为原生几何精度证明。

若只运行某个case的普通单项目输出，手动设置 `runFourPatternComparison=false` 即可。切换 `projectName="F192_AnodeBM_PlaneRGB"` 可恢复真实Anode/BM单元/阵列对比；其它旧项目不触发四图样分支。

### 真实Anode/BM平面波局部单元与阵列对比

选择 `projectName="F192_AnodeBM_PlaneRGB"` 时，该项目默认开启 `runPlaneWaveComparison`，运行后返回，不再重复执行普通单项目计算。计算由私有 `benchPlaneWaveComparison` 调用现有公共 `runCOEReflection` 完成。切换其它项目时此开关默认关闭；也可手动设为false，恢复当前项目的普通输出流程。采样基准开启时优先执行采样基准。

本项目保留现有Anode/BM图片的原始像素，只计算Anode反射与BM遮挡，不读取CF/IJP。原图1000×1000为4×4 Pitch完整输入单元，Pitch=59.06 μm；本次不裁切、更不推断最小晶胞。局部单元为4×4 Pitch（236.24 μm边长），阵列按config中的 `cfg.periods=[16,16]` 显式平铺为4000×4000（944.96 μm边长），均保留250点/Pitch。平面波、零额外高度；接收距离100 m仅服务于理想远场角分布，不代表实拍相机距离。

两组使用相同的±6°角域、1201×1201输出与接收像素积分。先分别计算620/550/450 nm强度，再用D65三点采样权重合成RGB代理；不是完整连续光谱色度，也不是相机预览。RGB直接使用同一个 `plotComposeReflectionRgb` 调用取得共同曝光和映射，禁止两图各自归一化。默认clipQuantile=0.9985、alpha=60、gamma=0.72；单波长图采用两组共同550 nm峰值作为0 dB参考、显示−60至0 dB。中心±1.5°图仅裁切显示，不重新计算，也不提高分辨率。

输出目录 `Output/F192_AnodeBM_PlaneRGB/PlaneWaveComparison/`：

- `plane_wave_comparison.png`：每行一组，依次为550 nm、全角域RGB、中心放大RGB。
- `local_rgb.png` / `array_rgb.png`：共同曝光的纯RGB图，无坐标；上方对应正θy。
- `geometry_and_profile.png`：实际有效反射区域与两组共同参考的550 nm水平剖面。
- `plane_wave_summary.csv`：各波长总反射率、窗口捕获率、零级占比、Parseval误差及物面面积。
- `plane_wave_comparison.mat`：`report.local/array` 内保存原始spectralEnergy、channelSignal、cfg、角坐标、入射面积/功率、输入周期和重复数，另保存共同显示映射、局部有效复反射图和汇总表；不保存无关大网格。

图中的强度为**每单位总入射功率对应的角像素能量**，适合比较能量分布形状；两组物面平均辐照度相同，但16×16阵列的入射面积/总功率是4×4单元的16倍。不要把归一化图直接解释为相同辐照度下的绝对接收亮度。代码检查共同角网格、完整平铺前后总反射率相同、每组仅传播一次高度状态及逐波长Parseval误差。以上检查不是物面或光谱收敛性证明。实拍比较优先看斑点排列、间距与包络；颜色及绝对强度还依赖真实光谱、材料和成像条件。

2026-11-16原生250点/Pitch实算：两组总反射率均为0.1575，最大Parseval绝对误差2.78×10⁻¹⁷；±6°窗口捕获局部单元反射能量的93.12%/93.87%/94.98%，阵列的93.93%/94.63%/95.61%（依次620/550/450 nm），窗口外能量未重新归一化。局部单元呈现较宽的彩色干涉光点，阵列呈现更窄的规则点阵，两组仍有明显横纵方向增强；与实拍有定性结构相似，尚未匹配其近圆形包络、光点宽度及蓝绿色颜色分布。三波长不能用来排除连续光谱影响。验证：新比较分支独立复制模块运行测试、原主脚本非零高度回归和D65色度/角域回归均通过（共3项）；新测试使用临时低分辨率夹具，仅验证流程/归一化/私有函数解析，不用于真实版图精度证明。

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

不再选择 simple/multilayer。`stack.layerNames` 同时决定读取哪些图片/材料；未声明的图层不读取。名称必须是唯一的有效 MATLAB 字段名，可自定义。`stack.layerLevels` 与列表一一对应，正整数按从下到上非递减排列；省略时按列表顺序逐项叠放。相同层级表示同一物理层中的互斥材料区域，其图片不允许交叠。旧Case示例的CF_R/G/B均为层级4，表示同层互斥。主脚本默认和F192/COEB3采用下述独立膜层配置，支持真实交叠。

实际RGB CF默认制作顺序为从下到上 **R→B→G**。直接重排名称列表即可自定义，不新增数字颜色编号；层级始终与列表一一对应并按从下到上递增：

```text
stack.layerNames = ["Anode", "IJP", "BM", "CF_R", "CF_B", "CF_G"]
stack.layerLevels = [1, 2, 3, 4, 5, 6]
```

例如改为R→G→B时，将最后三项改为 `"CF_R", "CF_G", "CF_B"`，层级仍为 `[1,2,3,4,5,6]`。材料参数按名称绑定，重排列表不改变该颜色的n/k或厚度。仅在掩膜实际重合处叠加膜层；缺失的膜不会保留空气占位层。交叠处各膜保持其配置厚度，局部总高度为所有存在膜厚之和，统一计算多界面反射、吸收及相位，并补足到共同参考面。额外`heightMap=0`不删除这些膜层相位。当前不模拟交叠处减薄或平坦化。

F192和COEB3均为1000×1000像素的4×4 Pitch单元，每Pitch 250像素；配置目标64×64 Pitch时各轴重复16次，得到16000×16000网格。按阈值0.5、BM反相后的原图统计，CF交叠分别占整个单元3.2416%和3.5168%，均位于BM实体上，且无三色共同交叠；BM上的交叠仍参与完整膜堆计算。

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

平面波可显式设置物面照明范围（长度单位m）：

当前16个项目config统一采用`type="plane"`、`integrationMethod="direct"`、`wavefront="plane"`，默认1×1 mm矩形照明。这里采用不同波长非相干叠加、单个入射模式内部空间相干的单方向近似，不称为“完全空间非相干平面波”。原ConvolutionRGB/ParaxialRGB项目名称及旧结果用于历史追溯；这些项目的当前config也已切换为plane，旧point的distance/diameter在plane模式下不生效。若要恢复有限发光面研究，应显式改回point和相应传播假设。

```text
cfg.source.type = "plane"
cfg.source.apertureShape = "rectangle"
cfg.source.rectangleSize = [0.001, 0.001]
cfg.source.circleDiameter = 0.001
cfg.source.centerXY = [0, 0]
```

`rectangleSize`是矩形[x全宽,y全高]；切换`apertureShape="circle"`时使用`circleDiameter`全直径。`none`为历史全物面均匀照明。该设置在物面上形成均匀强度足迹，不修改任何反射图片；每个波长内部仍空间相干，不等于完全空间非相干照明。完全空间非相干场不能同时限定为单一传播方向，需要另行规定角分布/互强度模型。

照明足迹与输入物面、`pupilRadius`共同取交集。平面波的incidentIrradiance为实际照亮交集内的均匀辐照度，incidentArea及入射功率分母使用该交集的离散面积；不会把变小的照明区重新放大到原物面总功率。`result.illuminationAperture`保存请求面积、实际面积、覆盖比例及是否超出物面。照明超出物面时打印提示，不自动扩大或平铺；若需完整1 mm光斑作用于大屏，应显式提供足够大的反射物面。例如case01~04的局部物面仅0.11812 mm宽，8×8 Pitch阵列也仅0.47248 mm宽；1 mm矩形/圆形照明均覆盖整个当前物面，切换形状不会改变这些小物面的结果。物面足够大、能看见照明边界时，矩形/圆形才会体现不同的衍射响应。

main当前默认`projectName="case01"`、`runFourPatternComparison=false`，直接执行所选项目的完整指标和绘图流程。四图样批量对比需手动设true。周期物面raster基准依赖全物面均匀照明，须显式设置`cfg.source.apertureShape="none"`。

2026-11-16已实际跑通case01的完整普通主流程，输出在`Output/case01/`；照明交集面积约0.223237 mm²。新矩形/圆形照明、入射功率、偏心光斑、空交集和超出物面不自动平铺等检查通过；完整仓库89项测试通过、0失败、0未完成，记录在本地`Output/IlluminationValidation/`。此次未更换夫琅禾费传播积分器，变更的是物面照明窗及对应入射功率归一化。

2026-11-16随后按用户要求，case01~04逐一关闭`runFourPatternComparison`、平面波、完整普通main流程重跑（均为8×8 Pitch），与旧`Output/FourPatterns/four_patterns_results.mat`中的array结果对照：四组逐波长spectralEnergy、总反射率、窗口捕获比例及零级比例的最大绝对差全部为0；使用旧共同曝光映射重绘后的RGB最大绝对差也全部为0。本次默认普通main对每个case单独确定曝光，因此原生RGB数组与旧共同曝光图的最大显示差分别约0.1450、0.1473、0.00429、0.00759（0~1范围），不代表物理强度改变。1×1 mm照明完整覆盖当前0.47248×0.47248 mm物面，覆盖率均为1。数值表为`Output/SingleFlowComparison/comparison_summary.csv`与`spectral_comparison.csv`，三行对比图为同目录`ordinary_vs_batch_rgb.png`，依次是旧共同曝光、本次独立曝光、本次按旧共同曝光重绘。绘图统一正θy向上；旧批量纯PNG使用flipud保存、普通main纯PNG直接保存，直接对照无坐标PNG时还须注意上下方向约定。原始旧批量结果未覆盖；运行结束后main恢复case01、四图样开关false。

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

### IJP：图片颜色由invert解释

实际项目默认顺序是 **Anode → IJP → BM → CF_R → CF_B → CF_G**；Case1~3保留显式同层互斥CF设置。IJP全平面连续覆盖，上下表面平坦，n=1.5、k=0、厚10 μm。Case1~3已提供与其它图层严格同尺寸、同原点的全白IJP.png/.bmp/.jpg，不再需要特殊的无图片输入流程。白色表示“此处有膜”，透明性由k=0定义；透明膜仍有界面反射和传播干涉，绝不能以删除复反射相位的方式模拟透明。

```text
stack.layers.IJP.file = "IJP"
stack.layers.IJP.method = "nk"
stack.layers.IJP.n = 1.5
stack.layers.IJP.k = 0
stack.layers.IJP.thickness = 10e-6
```

F192/COEB3的IJP原图为全黑，必须保持 `stack.layers.IJP.invert = true`，反相后为整面有膜；旧Case全白图采用false。invert只解释几何覆盖，不改变材料n/k、膜厚或光学相位。

IJP下方仅Anode图样范围内为反射底，其余位置沿用理想匹配吸收终端；覆盖IJP/BM/CF后的界面反射仍计算。工程接口保留coverage="full"便于算法验证，项目统一读取图片，并通过invert声明黑白含义。

厚膜干涉对光谱步长更敏感：仅10 μm、n=1.5的IJP在380 nm附近，往返相位变化周期约4.8 nm，5 nm可能不能分辨局部干涉条纹。当前常规方案比较按用户选择暂用统一5 nm；对很接近的候选排序或精细光谱形状，再用更细步长复核。主脚本保留厚膜采样提醒；完整膜堆的最严要求可能比单独IJP更高，建议步长不是精度保证。

<a id="convergence-decisions"></a>

### 收敛性结论、当前选型与验证记录

以下汇总用于后续参数选择与交接。**当前采用的计算参数与已经通过的收敛判据应分别理解**：64 Pitch、5 nm、±6°/1201、有限源5×5是目前结合计算资源和用户选择的常规配置，不代表所有物理量均已联合收敛。详细扫描条件、数据表和公式继续保留在本节后的采样/Pitch专项及输入单元说明中。

| 因素 | 已有验证结果 | 当前采用及适用范围 |
|---|---|---|
| 光谱步长 | 10 μm、n=1.5的IJP及完整膜堆均有快速干涉；完整膜堆5 nm的区域R曲线最大重建差约30.05个百分点，0.25 nm约0.37个百分点。宽带积分可比局部曲线更早稳定 | 方案常规比较统一5 nm；候选差异很小、材料/光源存在尖锐光谱或关注干涉细节时，再以更细步长复核。局部曲线误差不能直接当成方案差值或排序误差 |
| 输出角网格 | ±6°/801相对1601的共同角图L1差约占窗口能量9.04–25.70%，所选积分指标却可相差仅约0.0017个百分点；真实±1°/801比宽图裁切有更细采样 | 采用±6°/1201，间隔0.01°；属于当前计算预算的选择，不能据此宣称窄峰形状已充分分辨。积分量、角图和峰宽分别检查 |
| 像素内积分 | 在±1°/401、16 Pitch的已测条件下，最少3点相对7点的共同角图差约0.044–0.154% | 保留receiverBinSamples=3；它是每轴最少点数，实际值可能由生产算法自动提高，不等于输出角网格数 |
| 每轴Pitch数量 | 对同一无限周期效率参照，128 Pitch通过D65辐射/明视觉判据；192 Pitch通过所测81波长全部ROI判据，256 Pitch进一步确认；64 Pitch未通过该判据 | Case1/2/3暂用64×64 Pitch作为资源折中，原生64像素/Pitch时为4096×4096物面；需要上述宽带/逐波长效率标准时，分别用128/192 Pitch并复核资源 |
| 物面像素数 | 每Pitch 32点相对原生64点的共同角图差约1.99–3.21%；64→128点的硬像素保持细化差不超过约0.108%，后者仅验证同一栅格的积分 | 推荐64像素/Pitch起步；50 μm Pitch对应0.78125 μm间距。细边界验证用真实重新栅格化的128像素/Pitch输入，间距0.390625 μm；插值放大旧图不能证明几何收敛 |
| 有限光源采样 | 在550 nm、8×8 Pitch、距离0.1 m/直径0.2 mm圆盘条件下，3×3相对25×25参照的共同角图差约6.40%；17×17约0.070%；5×5、9×9误差不单调 | 按用户选择保留5×5（圆盘内实际21点）。当前样例照明仍为plane，该参数不触发源点循环；有限距离光源的结果不能直接套用平面波Pitch结论 |

效率的初始判据为`绝对差 ≤ 1e-4 + 0.01×abs(参照量)`，反射率按0–1计，1e-4即0.01个百分点，并非纯粹的“相对误差小于1%”。弱阶次使用同一绝对容差。最细样本的自比较为零不构成独立精度证明；推荐Pitch档位还要求更大样本继续通过。早期圆形ROI与Pitch专项的方形ROI需各自使用本阶段固定定义，不能逐行混比；自动S编号不能代替固定理论阶次。

自动平铺的正确性已于**2026-11-16**验证：42项COE回归通过，包括公共入口下2×2 Pitch、4×4 Pitch及2×4 Pitch输入单元与预先平铺完整图片的传播结果对照，以及目标不可整除、图层原始尺寸不一致的报错检查。Case1、Case2、Case3的实际config也逐项检查通过：输入16×16 Pitch、横纵各平铺4次、目标64×64 Pitch。此项证明图像装配和调用行为正确，**不替代物理采样收敛验证**。

记录阶段与复核入口：

- **2026-11-13五类采样扫描**：68组空间计算；基于当时16 Pitch Case1输入。原始数据在`Output/Case1/SamplingConvergence/`的`spectral.csv`、`ijp_spectral.csv`、`spatial.csv`、`spectral_reference.mat`、`sampling_convergence.mat`；保留半格偏移、材料/区域表、实际积分点数及源码快照。
- **Pitch专项**：数据在`Output/Case1/PitchConvergence/`的`pitch_summary.csv`、`pitch_convergence.mat`、`production_validation.csv`；保留无限周期参照、逐波长/宽带分开的通过标志、ROI及6组生产传播交叉检查。使用主脚本`runSamplingBenchmark=true`、`samplingBenchmarkMode="periods"`复核；五类扫描用mode="all"。
- **64 Pitch完整仿真**：数据在`Output/Case1_Pitch64/`，81个波长均已完成；`run_validation.json`记录原始数组/CSV、角域积分和完整输入哈希复核。此结果只代表实际运行的64 Pitch；准备好的128 Pitch图片项目未被标记为已完成完整仿真。
- **2026-11-16平铺回归与真实配置检查**：本地日志为`Output/repeat_config_check_20261116/matlab_verified.log`，回归记录为同目录`regression_results.mat`；对应测试在`tests/COEMultilayerTest.m`及现有COE回归套件中。

`Output/`默认不进入Git。复核或迁移时应单独备份所需MAT/CSV/日志及对应配置、原始输入和源码版本；本README中的已测条件与数值为版本化交接记录。2026-11-13的16 Pitch基准与2026-11-16的64 Pitch样例配置属于不同阶段，不能把旧输出当作当前config重新运行的结果。更换图样、材料、厚度、光源、目标角域或ROI后，应重新检查相关采样维度；离散Parseval通过仅说明离散能量核对，不保证光谱、空间或ROI都已收敛。

### 采样收敛基准

选定项目后，将主脚本 `runSamplingBenchmark=true`，完整运行。默认关闭；开启后复用项目配置、用户图片与统一膜堆，完成基准后返回，不执行常规绘图/导出。实现位于 `src/private/benchSampling.m`，不增加公共接口。基准要求图片输入、零额外高度、空外部复系数及无限光瞳；不调用archive，不改写Input。

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

“有限源采样”指照明光源发光面的积分，不是版图重复数。例如每轴17点是在发光圆盘的包围方形内建立17×17候选点，本次圆盘内实际225点；每点分别照明完整版图，再叠加强度。版图重复数由cfg.periods独立控制，历史有限源扫描为[8,8]、当时Case1为[16,16]；当前Case1/2/3的目标配置为[64,64]。

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

这是当前平面波、图样、材料及所列ROI的结论，不是所有光源/观察条件的通用Pitch数。有限距离点源的入射曲率和照明范围会随物面尺寸改变；选择5×5光源采样不消除这种物理变化。正式改用128/192 Pitch时，应准备对应的8192×8192/12288×12288完整同坐标用户图片，且核对远场条件与资源。原Case1仍保留16-Pitch输入，另有显式平铺的64/128-Pitch项目，见下节。

`pitch_summary.csv`保留逐波长和宽带的分开通过标志及推荐；`pitch_convergence.mat`保留各波长/尺寸/ROI原始效率、无限阵列参照、完整配置与源码快照；`production_validation.csv`保存交叉验证；`pitch_convergence.png`给出误差曲线，`pitch_images_550nm.png`在共同dB标尺下展示不同Pitch数的单色角图。当前交叉验证最大角图L1相对差约0.095%（生产代码采用有限midpoint积分），精确周期积分的最大Parseval差约7.5e-15。

### 完整图片项目与64 Pitch全流程结果

按用户要求，已基于Case1六张PNG原图显式准备两个完整项目，未缩放、插值或修改膜堆；各图层使用相同平铺方式和原点，采样间距仍为0.78125 μm：

| 项目 | 原图平铺 | 每层完整图片 | 完整物面 | 接收距离 |
|---|---|---|---|---|
| Case1_Pitch64 | 4×4 | 4096×4096 | 3.2×3.2 mm | 100 m |
| Case1_Pitch128 | 8×8 | 8192×8192 | 6.4×6.4 mm | 250 m |

各项目保留Case1的Materials、统一Anode/IJP/BM/CF参数和全白10 μm IJP。`input_manifest.csv`记录源文件/输出文件SHA-256、平铺倍数及覆盖率；逐像素验证完整输出与原图平铺一致。100/250 m分别超过最短波长380 nm下约54/216 m的常用夫琅禾费尺度；距离改变接收面尺寸和辐照度，不改变角域归一化反射效率。

用户随后允许适当降低Pitch数，本轮受内存资源限制，实际完整运行的是**Case1_Pitch64**；128项目已准备并保留，未将其声明为已完整仿真。主脚本当前选择`projectName="Case1_Pitch64"`、`useProjectConfig=true`、`runSamplingBenchmark=false`，可直接重跑。沿用用户确认的**plane照明**，`cfg.source.samplesPerAxis=5`保留为有限源参数，在本次plane模式下不启用5×5源点循环。

此次使用现有完整物面传播链，D65 380–780 nm、5 nm步长的81个波长全部计算，±6°/1201输出、额外高度0，仅一次传播状态；未使用专项基准的周期积分加速替代本次仿真。大尺寸多波长计算会打印逐波长进度。模型及数值算法未作变更。

| 指标 | 结果 |
|---|---:|
| 总反射率，辐射度 | 9.72003% |
| 总反射率，明视觉 | 9.59359% |
| 中心2.5°半角内反射率，辐射度 | 7.87733% |
| 中心2.5°半角内反射率，明视觉 | 7.67161% |
| ±6°窗口内辐射反射率 | 9.13206% |
| ±6°窗口捕获的总反射能量比例 | 93.95096% |

结果位于`Output/Case1_Pitch64/`，包括D65全范围/局部图、逐波长图、几何、累计曲线、接收面图、CSV及原始MAT。`coe_reflection_results.mat`约4.19 GB，含81×1201×1201的原始光谱能量及其它保存字段；`run_validation.json`记录独立核对。已按空间块读取全部81波长，检查非负有限能量、逐波长能量和、宽带圆锥积分和两套完整输入的像素/哈希一致性；原始数组与逐波长CSV最大绝对差4.86e-16，最大离散Parseval差1.78e-15。控制工具曾在等待期间超时，MATLAB在后台继续完成了全部导出；通过输出文件及原始数组确认完成状态，不能仅以工具超时推断计算失败。

本轮64 Pitch用于完整流程及图像观察，不改变上一节128/192 Pitch效率收敛结论。色彩图采用共享曝光和亮度增强，物理指标仍来自原始能量；当前是理想角分布，不是相机图像。

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

所有层必须覆盖同一个输入矩形、采用相同原点与像素网格，行为 y、列为 x。关闭平铺时，该输入矩形就是完整计算物面。图片本身不提供物理Pitch。默认仍按完整物面读取；用户可显式开启stack.geometry.repeatToFit，将所有图片作为同坐标输入单元平铺到cfg.periods指定的目标，默认不重采样、裁切或改写输入；仅显式指定inputSamplesPerPeriod且计算采样率更低时，允许在内存中进行共同像素中心最近邻降采样。

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

### 实际项目：物面采样率与Pitch数分开选择

F192/COEB3的原图采样率固定记录为250；计算采样率可以降低，但须先运行收敛基准。仅改变cfg.samplesPerPeriod而不声明原图采样率会改变物理尺度解释，因此应同时配置：

```text
stack.geometry.inputSamplesPerPeriod = 250
cfg.samplesPerPeriod = 125
```

这里125只是语法示例，不是已通过的推荐值。原图保持1000×1000，程序先在所有图层相同像素中心位置最近邻采样为500×500，再按原来的4×4 Pitch物理范围平铺。原始图片不会覆盖或裁切；不允许用该功能上采样。inputSamplesPerPeriod=0（默认）表示沿用cfg.samplesPerPeriod，不重采样。结果geometry记录inputGridSize、sampledInputGridSize、inputSamplesPerPeriod、samplesPerPeriod和repeatCounts。

主脚本设 `runSamplingBenchmark=true`、`samplingBenchmarkMode="raster"`，逐项目运行。基准直接读取完整4Pitch原图单元，扫描40/50/64/80/100/125/160/200/250像素/Pitch及0.25/0.5/0.75像素采样偏移，不先展开64Pitch物面。所有候选与原生250比较：380:5:780nm的总反射率、±6°方窗、2.5°圆锥及0.3°外非零级能量，以及同一组物理衍射阶次的强度L1误差。积分量容差为1e-4+1%×参考值，角分布L1容差1%；宽带D65辐射和视见函数积分另列。要求候选和更细候选的所有采样偏移均通过，避免幸运对齐。

若筛出较低采样率，再经完整4Pitch有限物面生产传播逐波长复核；这只能说明当前离散原图相对误差，不能证明连续几何已收敛，也不能代替64Pitch窄峰形状验证。若没有较低候选通过，则保留250，生产复核只检查550nm参考状态。结果在 `Output/项目/RasterConvergence/`，包括raster_summary.csv、production_validation.csv、MAT记录和误差/单缓冲区内存图。内存列是一个复数网格或regionMap的理论大小，不是总进程峰值。

再设 `samplingBenchmarkMode="periods"` 扫描每轴Pitch数。现在支持完整4Pitch输入单元，不再限定2Pitch晶胞；只考虑原图单元的整数重复，向下扫描到一个完整输入单元，最大值为当前cfg.periods。复用自相关精确积分并核对生产传播和Nyquist窗口Parseval；固定ROI内若包含多个无限周期阶次，应将全部阶次求和。逐波长与D65宽带的通过条件分别报告，ROI效率与峰高/峰宽不同，不能据积分通过宣称较小阵列具有相同峰形。

250像素/Pitch、64Pitch为16000²网格，仅一个complex double网格约3.815GiB；相同采样率下32Pitch降至8000²，单网格约0.954GiB，16Pitch约0.238GiB。实际流程同时保留多个实数/复数网格，因此优先依据指标收敛减少Pitch数，不能仅凭单数组大小推断整个任务可用内存。

### F192 / COEB3本轮采样结果

采用R→B→G、交叠区完整膜厚、黑色IJP反相，以及当前配置的占位材料常数。原64Pitch作为有限物面资源比较参照，不作为无限阵列或连续几何真值。380:5:780nm、D65，积分容差1e-4+1%×参考值，角分布强度L1容差1%。

| 项目 | 物面采样率结论 | 相对64Pitch的D65固定ROI最低通过周期数 | 当前共同配置 |
|---|---|---|---|
| F192 | 40～200未全部通过角分布门槛，保留原生250 | 32；最大归一化容差比0.7004 | 250像素/Pitch、32×32Pitch |
| COEB3 | 40～200未全部通过角分布门槛，保留原生250 | 24；32Pitch的最大归一化容差比0.5725 | 250像素/Pitch、32×32Pitch |

这里的固定ROI包括±6°方窗，以及所列阶次中心的0.08°/0.10°半宽方形ROI；D65辐射与视见函数积分均参与筛选。不是主流程2.5°圆锥指标或全角图逐像素精度的独立保证。125像素/Pitch的积分误差虽小，角分布仍未全部通过，因此未采用。

两套4～64Pitch都未通过全部逐波长ROI相对无限周期阵列的严格条件；最差处410nm附近的某些阶次紧贴固定ROI边界，有限峰宽使能量仍在边界内外分配。不得将资源推荐32Pitch称为全面收敛，也不得将64Pitch自比较作为精度证据。需要单色细节、峰形或新的实际n/k数据时应重跑相应基准。当前32Pitch仅用于D65固定ROI方案比较，图峰宽/峰高随物面尺寸变化。

结果位于各项目Output下的RasterConvergence和PitchConvergence，包含CSV、原始MAT、生产交叉验证及图像。Pitch表中原Recommended/BroadbandRecommended对应无限阵列；新增CurrentD65ToleranceRatio/CurrentD65Recommended仅对应本轮最大有限物面64Pitch。最大样本自比较不会被选为较低采样建议。

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

### 输入单元自动计算平铺次数

Case1、Case2、Case3现在统一使用目标64×64 Pitch、每Pitch 64像素、D65/5 nm、±6°/1201和100 m接收距离；保留plane照明和有限源每轴5点。64 Pitch是目前跑通完整流程的资源折中，不改变128 Pitch宽带效率收敛建议。原始PNG/BMP/JPG文件不修改。

```text
stack.geometry.repeatToFit = true
cfg.periods = [64, 64]
cfg.samplesPerPeriod = 64
```

已知每Pitch像素数后，输入单元的横纵Pitch数分别为`[图片列数,图片行数]/samplesPerPeriod`；重复次数为`cfg.periods./输入单元Pitch数`，输出像素数为`cfg.periods*samplesPerPeriod`。无需人工填写repeat次数，不从图案外观猜Pitch，也不搜索图片内部的最小数学周期。现有1024×1024图片作为16×16 Pitch单元，横纵各重复4次。

| 输入单元（横×纵） | 每Pitch像素 | 图片列×行 | 目标 | 重复次数（横×纵） |
|---|---:|---|---|---|
| 4×4 Pitch | 64 | 256×256 | 64×64 Pitch | 16×16 |
| 2×2 Pitch | 64 | 128×128 | 64×64 Pitch | 32×32 |
| 2×4 Pitch | 64 | 128×256 | 64×64 Pitch | 32×16 |

所有真实图层图片必须同尺寸、同原点，包括IJP。默认原图行列须为samplesPerPeriod的整数倍；显式降采样时按下述inputSamplesPerPeriod解释原图；目标Pitch数须可被输入Pitch数整除，不满足时报错。相同边界的材料、缺陷也会一起平铺，因此用户应提供能首尾衔接的完整光学晶胞，不含重复首尾端点。heightMap仍表示完整目标物面，不随图层自动平铺。MAT输入不启用此开关。

默认repeatToFit=false，保留原完整物面输入行为；已有Case1_Pitch64/128完整项目无需开启。运行时打印输入Pitch数和[x,y]重复次数，保存的stackGeometry或公共接口result.layerGeometry中有inputGridSize、inputPeriods、repeatCounts记录。该授权的显式平铺不引入图形生成器，不改变统一膜堆，也不增加相机模式。

2026-11-16验证：42项COE回归通过，包含2/4 Pitch、矩形输入、目标不可整除和图层尺寸不一致；Case1/2/3真实配置均确认输入16×16 Pitch、重复4×4、目标64×64 Pitch。

每Pitch默认推荐**64像素**：50 μm Pitch对应0.78125 μm间距，适合作为当前方案比较的起点；128像素/Pitch对应0.390625 μm，适合细边界和窄间隙复核，但应提供相应真实分辨率原图。把64像素原图插值放大到128并不证明几何收敛。以目标64×64 Pitch为例，两档完整物面分别为4096²和8192²点，128档的物面数组约增加到4倍。每Pitch像素数不是衍射角网格，二者应分别检查。

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
- R/G/B CF 覆盖各自对应的 BM 开口，其区域大于 BM 开口；实际项目允许各色CF交叠，并按制作顺序叠膜；同层互斥仅适用于显式采用相同layerLevels的示例。
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

各位置先按局部平行平面膜堆计算随波长变化的复反射系数，再进入指定频率网格的傅里叶积分传播 / 角域 / CIE 后处理。CF 局部厚度固定，位于 BM 实体上方时叠加在 BM 顶面，在开口内则从裸底平面起算。以各物理层级最大nk厚度之和为共同空气参考面（R→B→G独立层级配置为 `IJP厚度+BM厚度+R厚度+B厚度+G厚度`；同层互斥示例为 `IJP厚度+BM厚度+最大CF厚度`），补足各柱体顶面到该平面的往返相位。未包含侧壁或邻近柱体的矢量耦合，不能当作 RCWA/FDTD。

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

2026-11-16提交前执行完整仓库 `runtests('tests')`：85项通过、0失败、0未完成。覆盖四图样/真实Anode-BM平面波比较、面源卷积、图片平铺与采样、膜堆、色度、能量和独立模块运行。完整结果保存在本地 `Output/CommitValidation/validation_results.mat` 及CSV（输出不提交到Git）。该回归验证数值和程序行为，不等于实拍拟合或产品材料精度验证。

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

私有函数按用途使用三类前缀，文件列表按名称排序即可分组；分类不取决于是否被主脚本直接调用。主脚本 `COEReflection_main.m`、公共接口 `coeReflectionDefaults.m` 和 `runCOEReflection.m` 保留原名，文件内局部函数也保留原名。

| 前缀 | 文件数 | 用途与示例 |
|---|---|---|
| `core` | 18 | 正式仿真的输入、进度、物理计算、色度和指标，例如 `corePrepareLayerStack`、`coreReflectFilmStack`；`coreSynthesizeColor` 保留原始 XYZ 并生成色度显示 |
| `plot` | 5 | 绘图、样式、标记及 RGB 代理预览，例如 `plotFormatFarFieldAxes`、`plotComposeReflectionRgb` |
| `bench` | 6 | 可选采样及结构对比：`benchSampling`、`benchPeriods`、`benchRaster`、`benchInputCell`、`benchPlaneWaveComparison`、`benchFourPatterns`；可重复运行，不代表一次性工具 |

当前 `src/private/` 保留29个函数文件。只在一处使用的简短辅助逻辑直接并入调用处：反射默认值集中在 `coeReflectionDefaults.m`；膜堆默认设置和采样提示在主脚本；图片查找和材料光谱在 `corePrepareLayerStack.m`；两轴傅里叶积分在 `corePropagateReceiver.m`，两轴共用的一维变换保留为局部函数。按主脚本及两个公共接口的静态调用关系检查，没有完全不可达的独立函数；条件分支中的函数并非每次运行都会调用。

| 分类 | 函数 | 触发条件或作用 |
|---|---|---|
| 默认配置 | coreLayerDefaults、coreSourceDefaults | 图层模板和光源默认值被多处复用 |
| 运行进度 | coreCreateProgress | 主脚本 showProgress=true 时创建进度窗，通过回调接收逐波长/源点进度；不参与物理计算 |
| 文件读取与膜堆构建 | coreLoadProjectConfig、coreLoadSpectrum、corePrepareLayerStack | 按项目配置读取实际图样和材料；图层/光谱方式决定具体分支 |
| 照明与传播 | coreBuildIllumination、coreSimulateHeightStates、coreSimulateDiffraction、corePropagateReceiver、coreReflectFilmStack | 当前统一计算链；零高度只计算一个高度状态 |
| 指标 | coreComputeAngularMetrics、coreComputeDiffractionSpots | 角域、光斑和能量评估；采样提示直接在主脚本中给出 |
| 独立采样基准 | benchSampling、benchPeriods、benchRaster | 主脚本runSamplingBenchmark=true时运行，mode="all"隔离五类因素，mode="periods"评估每轴Pitch数，mode="raster"评估物面采样率 |
| 绘图样式 | plotFormatCumulativeAxes、plotFormatFarFieldAxes、plotMarkDiffractionSpots、plotStyleFigure | 当前图形显示，标记开关只改变显示内容 |
| 条件色彩处理 | coreSynthesizeColor、plotComposeReflectionRgb | 前者用于D65/custom/monochromatic；后者仅RGB代理模式 |
| 历史批量接口兼容 | coreCustomPattern、coreNormalIncidenceReflection | 主脚本不调用该旧几何/单界面分支；仅传cfg的历史runCOEReflection接口和兼容测试仍可能使用，当前不能算完全无用函数 |

以下5个文件已移至 `archive/input_generation/`，不再放在当前private目录：

| 归档文件 | 原用途 | 当前状态 |
|---|---|---|
| createCoeDiamondGeometry.m | 程序生成固定Diamond几何 | 已有用户图片，不在主流程调用 |
| createCoeImageExamples.m | 导出PNG/BMP/JPG示例输入 | 开发阶段数据准备工具，运行仿真不需要 |
| createCoeHeightExample.m | 生成可选二值高度条纹 | 演示工具；当前仅加载用户提供的heightMap |
| coeStackDefaults.m | 旧Anode/BM/CF固定结构默认值 | 已由主脚本中的统一图层默认设置替代，归档前仅被旧工具和旧测试使用 |
| prepareCoeStack.m | 旧固定COE图样与膜堆入口 | 已由 `corePrepareLayerStack` 替代，旧测试已迁移至统一图层接口 |

归档保存原始实现以便查阅，不加入MATLAB路径、不由主脚本/测试/打包入口调用。归档函数保留其原依赖，不承诺单独可运行；如需恢复，应在隔离副本中结合Git历史恢复相应上下文，不能直接将归档目录加入当前工程搜索路径。当前`src/`与`Input/`即可独立运行，不需要archive目录。

测试中的小尺寸图片和高度矩阵属于临时校验夹具，不是面向用户的输入生成程序；它们只在测试临时目录中创建，目的是验证实际文件读取、角域积分与回归行为，不覆盖用户输入。

### 程序责任与后续维护

| 文件或目录 | 用途 |
|---|---|
| src/COEReflection_main.m | 输入、配置、两状态计算、展示和保存 |
| src/coeReflectionDefaults.m、src/runCOEReflection.m | 保留兼容的公共配置/批量计算接口 |
| src/private/coreLoadProjectConfig.m | 项目文本配置解析 |
| src/private/corePrepareLayerStack.m、coreReflectFilmStack.m | 图样校验、区域查找表和相干膜堆 |
| src/private/coreSimulateDiffraction.m | 复振幅、接收面积分与基础指标；独立FFT能量诊断 |
| src/private/corePropagateReceiver.m | 二次相位卷积傅里叶积分、传播系数、像素内强度积分 |
| src/private/coreLoadSpectrum.m、coreSynthesizeColor.m | D65、CIE积分和色度显示 |
| src/private/coreComputeAngularMetrics.m、coreComputeDiffractionSpots.m | 角域与光斑积分 |
| src/private/plot*.m | RGB 代理预览、图形样式和辅助标记 |
| Input/Colorimetry/ | 原始CIE数据、元数据、许可与校验值 |
| AGENTS.md | 本模块维护规则 |
| dist/ | 历史打包工具与旧快照，继续工作以当前源码为准 |
| archive/input_generation/ | 已退役的几何/图片/高度生成和固定结构入口，仅作历史源码存档 |

src 顶层仅保留主脚本与两个兼容公共接口。仅在一处使用的简短逻辑优先直接写在调用处，不为减少单个文件行数而新增函数；有复用价值或需要独立验证的算法再保留函数边界。新增独立内部函数直接放在 private，使用小驼峰命名，文件名与首函数名一致；用户主流程和必要中文注释留在主脚本。现有数值算法回归由 tests/private/addCoeAlgorithmFixture.m 在临时副本运行，生产目录的私有解析由独立模块测试覆盖。

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
- 同层交叠报错：核对layerLevels及invert，相同层级表示互斥区域，不能含交叠。对于实际R→B→G叠膜，按本说明设置独立层级并保留掩膜和各层厚度；不要裁掉交叠区。
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


### F192实拍、旧孔径示例与当前D65图的差异诊断（2026-11-16）

用户确认：实拍为F192，手机照明光源有效发光面约5mm，到屏幕距离约5cm；可视为正对屏幕，目视效果接近照片。当前F192配置仍为plane、D65/5nm、32Pitch，不能把这张理想角分布直接当作该实验的预测。

本轮只读已有完整F192 MAT结果，没有启动新的传播。消融脚本及图在 `Output/F192_ImageDiagnosis/`：

- `spectral_display_ablation.png`：同一复反射场，A为81波长D65+XYZ+现有对数映射；B为同一D65结果关闭对数映射；C只取620/550/450nm并逐通道归一化，复现旧式RGB显示逻辑；D为550nm单色。四图统一显示±3°，实拍没有角度标定，不能按像素直接配准。C不是D65定量结果，D不是白光结果。
- `coherent_mean_diagnostic.csv/.png`：对实际区域复系数计算无限周期零级效率。550nm总反射率7.9104%，零级约0.5003%，占反射功率6.324%；人为去掉空间相位差会将零级改成7.3971%。后者仅用于说明相位敏感性，不能作为正式物理解法。实际零级弱需要检查真实n/k、膜厚、连续界面及其相干关系，不能给相位强制清零或人为加中心白点。
- `region_phase_550.csv`和`diagnostic_summary.json`保留计算依据。旧孔径示例为实数振幅掩膜、三个波长分别归一化和截断；当前为各局部膜堆的复反射系数及完整光谱XYZ积分，两者不等价。旧脚本的输出二次相位与当前省略该项的差异不会改变单独观察面的强度，不能靠补这项消除彩色射线。

光栅的同一非零阶角度随波长变化，因此连续谱叠加会形成径向色散轨迹，三条谱线则主要留下三组离散图样。当前toneStrength=3000、luminanceGamma=0.75又把弱尾部大幅提亮。减弱显示增强会隐藏弱尾部，但不改变实际能量分布。

5mm有效面源在5cm距离的全张角为2*atan(2.5/50)≈5.725°；550nm/59.06μm对应一个Pitch的角间隔约0.534°（完整光学晶胞可能有更细的阶次）。同时1.88992mm物面从中心到边中点的入射球面波光程差约8.929μm，即550nm下16.23个波长，平面波近似不适用。面源按互不相干源点的强度求和，空间亮度分布仍需实验约束。若均匀5mm面源把点阵洗平，这属于模型与实拍条件待校核的结果，不能擅自缩小有效发光面来制造相似图。

代码调整顺序：

1. `coreBuildIllumination.m`已有有限距离point源，输入应记录source.type="point"、distance=0.05、diameter=0.005；plane下distance与diameter不参与该计算。先做小成本单源及面源收敛诊断，不直接运行8000²、81波长、大量源点的暴力组合。默认每轴5点只用于演示，不能证明5mm面源已收敛。需要加速时，可在已验证近轴及移不变条件下，将带球面曲率的单源结果与角向辐亮度分布做逐波长强度积分/卷积，并与显式源点结果核对；不要对RGB增强图做模糊来替代物理源积分。
2. `coreLoadSpectrum.m`已有custom SPD接口；实拍应使用实际手机闪光灯SPD，而不是默认D65或人为三条谱线。目视与照片接近意味着不能把差异主要推给相机白平衡；照片本身仍不足以反演真实SPD。
3. 核对接收光路。`corePropagateReceiver.m`当前输出是远场角分布，cfg.propagationDistance只给出角度到接收面坐标与辐照度标定，改成0.05并不会切换到Fresnel近场或眼/镜头成像。只有观察几何确定后，才能决定是否另加有限距离传播或成像模型；“正对”本身不等于远场。
4. 实际材料CSV/膜厚、IJP所在位置和连续界面需核实。当前nkFile为空，所有n/k均为占位常数；不能把材料名CF_R/G/B等同于实际RGB透过光谱。保留R→B→G交叠及各层膜厚，不因视觉相似度改变几何或去掉相位。
5. `coreSynthesizeColor.m`继续以原始XYZ作为物理结果，同时查看线性曝光与温和显示映射；旧式RGB三线图只作排查。统一角度范围、曝光与标注后比较，不把显示效果相似当作物理验证。

本轮没有修改正式物理配置或传播算法。下一步应先匹配照明与接收条件，再用实测SPD和材料数据校核；完整大网格计算保持停止。


### Anode/BM三波长与非相干面源卷积基线

图样对齐入口为 `projectName="F192_AnodeBM_ConvolutionRGB"`，对照项目为 `F192_AnodeBM_ParaxialRGB`。两者复制F192的原始Anode/BM图片，维持250像素/Pitch和32×32Pitch；只声明Anode、BM，不读取CF/IJP。Anode采用统一反射率1，BM先作为反射率0的理想不透明黑区，厚度0表示有效边界。此处15.75%的总反射率来自阳极面积填充率，不能当作真实器件反射率。原F192/COEB3的膜堆和材料配置不修改。

三波长为620、550、450nm，RGB代理相对权重由共享CIE D65表在三处采样并归一化，约为[0.28404,0.33699,0.37897]。这是图样代理，不是完整D65色度积分；后续仍需恢复实际光谱与材料。

```text
cfg.source.type = "point"
cfg.source.distance = 0.05
cfg.source.diameter = 0.005
cfg.source.integrationMethod = "angularConvolution"
cfg.source.wavefront = "plane"
```

`integrationMethod="direct"`、`wavefront="spherical"`是原有默认：精确球面距离与辐射包络，各源点按强度求和。新增波前选项仅作用于point：

- `plane`：各点产生倾斜平面波，中心源对应光栅固有远场PSF；可表达在角度坐标中光源形状与光栅强度PSF的卷积模型。
- `paraxial`：保留中心源的二次相位exp(i*pi*r²/(lambda*d))，同时将偏轴源作用表示为倾斜相位。源点物面振幅视为相同，不使用精确距离衰减和beamHalfAngleDeg包络。

`angularConvolution`只允许point及显式plane/paraxial波前，拒绝将精确球面源静默替换为近轴源。均匀圆盘光源在每个波长按相同振幅、互不相干、平移不变假设进行强度线性卷积；核角半径取D/(2d)。5mm/5cm对应约2.865°半径，而不是5mm直接当接收面像素长度。source.samplesPerAxis在该快速方法中不控制源积分，核网格由输出角度步长决定，仍需检查输出角网格的收敛。

实现先将单源传播计算角域向外扩展至少一个源半径，再用归一化圆盘核线性卷积，最后裁切原观察窗。这样保留窗外单源图样卷入窗口的贡献，且没有FFT环绕；只归一化光源概率核，不把裁切后的反射能量补回1。补边后的角域不得超过近轴10°限制。每个波长保留result.sourceConvolution中的pointEnergy、pointClearEnergy、kernel、kernelStepRad、paddingBins、massError和assumption；原有spectralEnergy是卷积后的物理能量数组。

这两个选项是不同的光路近似，不是自动完成了眼睛或镜头的有限距离成像。光源5mm有效尺寸本身不能确定其空间角辐亮度、封装光学和空间相干性。若均匀非相干圆盘洗平点阵，应先核对源模型与实际光路，不能擅自缩小有效发光面来制造相似结果。

验证包括相同近轴假设下的显式源点求和、偏心亮峰附近的观察窗补边、零直径退化及核能量守恒。显式源点与卷积在所选对齐网格上的能量数组绝对差小于1e-11；另通过仅复制模块的公共入口验证真实private解析。此数值一致性不证明精确球面光路与近轴模型等价。

本轮结果写入各项目Output，统一图样比较写入 `Output/F192_AnodeBM_PatternComparison/`。比较图同时保留各图独立曝光（仅看形状）与共同曝光（看相对强弱），颜色标为三波长代理。


本轮实算已完成：两组均为8000×8000物面、±6°/1201输出及三个波长，运行含保存分别174.36秒、166.56秒。光源核575×575、角半径2.86479°，单源传播计算先扩展到±8.88°再卷积裁切；卷积质量检查误差不超过1.11e-16。35项COEMultilayerTest及独立复制模块的卷积入口测试通过，相关私有函数静态分析无问题。

结果见[图样独立曝光对比](../Output/F192_AnodeBM_PatternComparison/patterns_shape_comparison.png)、[共同曝光对比](../Output/F192_AnodeBM_PatternComparison/patterns_shared_exposure.png)及[原始指标](../Output/F192_AnodeBM_PatternComparison/pattern_metrics.csv)。独立曝光图只比较形状，不比较绝对亮度。点源远场基线有清晰离散点阵；均匀5mm非相干圆盘卷积将点阵显著洗平，保留入射曲率后同样不能恢复实拍点阵。因此这组简化假设尚未实现实拍图样匹配，后续需约束光源的空间角辐亮度、封装光学及相干性，不能把5mm有效发光直径直接当成已知完整源模型。


### 均匀圆盘未复现实拍后的代码调整方向

均匀圆盘卷积结果不能直接否定强度卷积形式；它只能说明当前源模型与观察模型尚未匹配实验。此前数值测试验证的是相同假设下的数学一致性，不是手机闪光灯出射场的物理标定。

只读敏感性脚本 `Output/F192_SourceModelDiagnosis/source_core_halo_scan.py` 使用已保存的Anode/BM点源结果与5mm/5cm圆盘结果，测试有效角分布中“窄角高斯核心+宽角背景”的功率混合；核心sigma=0.03°及功率占比都是假设扫描值，没有拟合实拍，也没有写入正式source配置。各图独立曝光，只看图样，不能比较绝对亮度。窄核只作用于已有±6°数据，展示限制在±3°，留足边界缓冲；没有重新运行传播。

在自由传播、均匀圆盘、源面空间不相干的近轴假设下，van Cittert–Zernike模型给出mu(b)=2*J1(v)/v，v=pi*D*b/(lambda*d)。D=5mm、d=5cm、lambda=550nm时，首个相干零点间隔约6.71μm；在59.06μm间隔处|mu|≈0.00814。该结论属于这个源假设，不是对实际手机闪光灯相干度的测量。封装光学、位置—角度相关性、接收选择和非均匀辐亮度都需要实际约束。

推荐的实现顺序：

1. 将`corePropagateSourceConvolution`中硬编码圆盘核抽离为独立源核构建步骤，支持经标定的二维有效角辐亮度图、角度采样与中心位置。5mm保留为空间出光口径，不能在未知封装光学时把它直接当作完整角分布。窄核心/宽背景仅可作参数敏感性诊断，不能因图样好看就宣布标定完成。
2. `coreBuildIllumination`逐步支持由实际源和封装光学产生的复场/相干模式，逐模式传播后按功率非相干求和；保持内存中只处理一个模式。独立LED源点不能改成固定相位的相干振幅相加。若有效角核模型足以解释实验，应先保留简单模型，不引入任意相干长度拟合参数。
3. 在`corePropagateReceiver`旁增加独立有限距离观察支路，按需要包含屏幕到入口瞳的Fresnel传播、有限瞳孔、镜头相位及对焦，再传播到观察面。保留farField为基准。正对屏幕、目视与相机相近，都不能单独确定该支路参数；仅调整z或倍率不保证恢复被源平均抹掉的点阵。
4. 点阵间距、包络、中心及非零级能量、对比度分别验证。先维持Anode/BM和三个波长，源核与光路匹配后再恢复材料和颜色。可用相同光路下普通镜面的光源反射像约束有效源像/观察响应，并用光源距离变化对比度检查角展宽假设。

本轮未修改生产算法或正式物理参数，只保留上述诊断脚本、原始CSV和对照图。一个未标定的照片不能唯一确定源相干性、有效角核和镜头参数，不能保证任意一组参数可以物理可信地复现实拍。
