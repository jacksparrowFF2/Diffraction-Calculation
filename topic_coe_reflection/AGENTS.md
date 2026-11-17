# COE 模块维护约定

- 用户运行入口为 `src/COEReflection_main.m`，集中呈现输入、参数、计算、展示及保存流程，并保留必要中文注释。
- `src/` 顶层仅保留主脚本，以及已供外部代码调用的公共接口 `coeReflectionDefaults.m`、`runCOEReflection.m`。未经明确需求，不新增公共接口。
- 新增内部计算、输入读取、绘图与保存辅助函数时，直接创建在 `src/private/`，不要先放在主脚本同级再要求用户手动移动。
- 内部函数及其文件统一遵循文末“函数命名与拆分规则”，采用 `core`、`plot`、`bench` 三类前缀。
- 保留 MATLAB 特殊目录名 `private`，不要对其单独 `addpath`，不要使用 `genpath` 展开私有目录。模块只需添加 `src/`。
- 移动函数时检查基于 `mfilename('fullpath')` 的输入路径，并同步修改所有调用、测试、说明与打包脚本。通过主脚本和独立复制模块的回归验证真实私有函数解析。
- 新的功能测试优先通过公共入口验证。现有数值算法测试的临时副本由 `tests/private/addCoeAlgorithmFixture.m` 管理，不为测试暴露生产私有函数。
- 实际版图按完整物面导入，不隐式缩放、平铺或改变阳极反射范围；额外 `heightMap` 默认保持 0 nm。
- 用户图片集按 `Input/项目名称/` 管理，主脚本使用 `projectName` 切换；输出对应 `Output/项目名称/`。CIE 数据共用 `Input/Colorimetry/`。不要在日常运行时覆盖用户输入；生成工具归档于archive/input_generation，不在主流程提供生成开关。

- 正式仿真不提供像素轮廓调整，用户输入阳极、BM、RGB CF 图片即可；不可重新暴露圆角、倾角、菱形尺寸等参数。参考图中实线是 PDL（本模型阳极），虚线是 Mask，不是 BM。预生成工具仅在archive/input_generation保留历史快照，不加入运行路径。
- 常用参数按项目保存为 `Input/项目名称/config.txt`，在光谱和几何计算前加载；配置覆盖默认值。解析文本不得执行 MATLAB 表达式；未知/重复/类型错误应报错。结果中保留配置记录。

- 使用说明、输入/配置规则、模型与指标解释、迁移交接统一维护在 `docs/README.md`；不再新增平行功能说明。Input 中的简短项目清单和 CIE 数据署名/许可保留在原位。

- 主脚本统一按stack.layerNames/stack.layerLevels构建膜堆，不再提供simple/multilayer选择。每项method为nk有限膜或reflectance不透明有效边界；默认RGB CF按从下到上R→B→G独立层级叠放，制作顺序通过layerNames列表自定义，交叠区保留各层膜厚；同层级仍要求互斥，旧示例可显式使用同层互斥配置。无基底膜用medium出射介质，COE外部默认absorber终端，不能把悬空膜未反射能量直接称为吸收。

- heightMap全零时只运行一次传播、单次色度和角域积分，并使用单状态图表及导出；非零高度保留两状态比较。历史baseline等数组字段可复用作兼容，不得重新计算相同结果。膜厚本身的相位不可因关闭额外heightMap而删除。

- 照明支持plane、有限距离/尺寸point和近轴TEM00 laser；有限源点按强度叠加，勿相干叠加。光谱与空间模型独立，支持D65/自定义SPD/单色谱线，不增加相机预览。入射功率按物面平均辐照度标定。
- IJP由图片及invert共同定义覆盖；旧示例为全白且invert=false，F192/COEB3为全黑且invert=true，与其它图层同尺寸同坐标，n=1.5、k=0、厚10μm，在Anode/BM之间。透明不代表无界面反射或无相位。厚膜需报告光谱欠采样风险。常规方案比较暂用5 nm；默认角域±6°、输出1201，有限圆盘源每轴5点；只裁切显示不等于提高计算分辨率。

- 当前仿真仅读取用户图样与可选heightMap输入；Diamond/图片/高度演示生成器及旧固定COE入口已存档，不再调用。回归测试应使用统一图层公共接口和临时输入夹具，不依赖archive。`coreCustomPattern` / `coreNormalIncidenceReflection` 暂为历史批量接口兼容保留，不得误认为无引用后直接删除。

- 采样基准由主脚本runSamplingBenchmark显式开启，默认关闭，内部实现放private。复用项目图片和统一膜堆，不改写Input、不增加相机或公共接口。重复数扫描先验证完整晶胞周期，仅在基准内存中显式平铺；物面抽取与硬像素保持细化要分别标记，不能把插值图当作更高分辨率真实版图。光谱曲线误差、宽带积分、固定ROI和角分布误差分别记录，最细参照样本的自比较不作为独立精度证明。
- samplingBenchmarkMode="periods"专门评估每轴Pitch数：以完整晶胞自相关和有限阵列三角权重精确积分矩形ROI，必须通过生产传播交叉验证并检查Nyquist窗口Parseval。平面波下与无限周期阶次效率比较；D65宽带通过与逐波长通过分开记录。不得将平面波推荐Pitch数直接宣称为有限距离光源的通用结论，不因推荐数改变现有用户图片尺寸或隐式平铺正式输入。

## 函数命名与拆分规则

本规则适用于本 COE 模块中新增或重命名的独立私有函数。名称采用 lowerCamelCase（小驼峰），以用途前缀开头，后接明确的功能名称；文件名必须与首函数名一致。

| 前缀 | 用途 | 示例 |
|---|---|---|
| `core` | 正式仿真的输入读取、配置、膜堆、照明、传播、色度及指标计算 | `corePrepareLayerStack.m`、`coreReflectFilmStack.m`、`coreSynthesizeColor.m` |
| `plot` | 绘图、坐标轴格式、样式、标记及 RGB 代理预览 | `plotFormatFarFieldAxes.m`、`plotMarkDiffractionSpots.m`、`plotComposeReflectionRgb.m` |
| `bench` | 可选的采样、收敛性及精度对比验证 | `benchSampling.m`、`benchPeriods.m` |

- 按职责分类，不以是否被 main 直接调用为依据；间接参与正式仿真计算的函数同样归入 `core`。
- `bench` 表示可选验证，可以重复运行；真正的一次性输入生成工具继续归档到 `archive/input_generation/`。
- 主脚本 `COEReflection_main.m` 和公共接口 `coeReflectionDefaults.m`、`runCOEReflection.m` 保留原名及调用兼容。
- 文件内局部函数不增加分类前缀，跟随所属文件管理。
- 仅在一处使用的简短辅助逻辑优先直接并入调用处，不为缩短文件而单独创建函数；有复用价值或需要独立验证的算法保留函数边界。
- 重命名时同步修改文件名、首函数声明、调用点、测试和说明，并检查打包脚本；完成回归测试及独立复制模块运行验证。

- 用户已授权显式图片单元平铺：stack.geometry.repeatToFit默认false，Case1/2/3设true并以cfg.periods=[64,64]为目标。根据原图列/行数与samplesPerPeriod计算[x,y]输入Pitch和重复数，不猜最小周期；各图层原始尺寸一致、目标必须整除，不缩放/裁切/改写输入。额外heightMap仍覆盖完整目标。记录inputPeriods/repeatCounts，保持默认完整物面行为和MAT入口兼容。

- F192/COEB3原生图片为4×4 Pitch、每Pitch250点；允许显式stack.geometry.inputSamplesPerPeriod=250配合较低cfg.samplesPerPeriod，在内存中按共同像素中心最近邻降采样后平铺，不覆盖输入。默认0沿用计算采样率。先运行raster基准，再独立评估periods；不得用原生参照自比较宣称连续几何收敛。periods基准使用完整原图单元，ROI内无限阶次须全部积分；降低周期数的宽带/逐波长结论及峰形变化分别说明。

- 图样对齐可显式选择F192_AnodeBM_ConvolutionRGB/ParaxialRGB：仅Anode/BM，BM为理想不透明黑区，不读取CF/IJP；三波长按D65采样权重，仅称RGB代理。原F192/COEB3完整膜堆不改写。
- 面源卷积为source.integrationMethod="angularConvolution"，仅point源并显式选择wavefront="plane"或"paraxial"，假设均匀非相干圆盘及源点相同振幅。球面精确直接积分保留direct/spherical默认，不得静默替换。逐波长在物理强度上做线性卷积，先向外补足源角半径再裁切，核积分为1但不能重标定最终窗口能量。plane对应远场PSF平移模型，paraxial另含中心源二次相位；两者均不是完整相机/目视光路。验证应与相同波前假设下显式源点求和对照。
- plane照明范围由source.apertureShape=none/rectangle/circle、rectangleSize全宽高、circleDiameter全直径和centerXY定义，单位m，不改变反射图样。与物面/几何光瞳取交集，incidentArea及功率分母使用实际照亮交集；照明越界只报告，不自动扩大或平铺输入。API默认none保留兼容，各项目当前显式1×1mm矩形、plane/direct。单波长模式内部仍空间相干，不同波长/独立源模式叠加强度，不得称为完全空间非相干平面波。main默认case01且runFourPatternComparison=false；批量比较需显式开启。

- 用户2026-11-17最新统一要求覆盖旧项目运行参数：全部24个Input config采用4×4 Pitch、接收面25cm、±5°、空间相干plane/direct、无独立人眼入瞳，D65的380:10:780nm共41波长，CIE XYZ→sRGB合成。各项目Pitch及膜堆保留，计算采样64点/Pitch；旧config冻结在tests/fixtures/coe_legacy_configs用于历史回归。用户已明确全部4×4，故对大于4×4的图片可显式cropToFit选中心整周期；pixelReplication只在显式开启且整数倍时复制原始像素块，不能宣称提高真实轮廓分辨率。原图片不改写，记录原周期、选取范围、重复数与复制倍数。两开关默认false，禁止静默裁切或默认上采样。主脚本仍默认case01，不自动启用RGB专用历史基准。
