# UDC 透射衍射成像 V2

在 MATLAB 中运行 [UDCImagingV2_main.m](../src/UDCImagingV2_main.m)。默认读取 [case01 配置](../Input/case01/config.json)，以同一张透射掩膜和 V1 的 `USAF500.png` 场景运行 `specifiedTransmittance`（给定功率透过率）和 `nk`（按薄膜 n、k 与厚度计算）两种配置。主脚本的 `transmissionModes` 可只选择其中一种。配套运行 [COEReflection_main.m](../../topic_coe_reflection/src/COEReflection_main.m) 可得到同一版图的反射衍射。

- `Input/<项目>/` 存放透射掩膜和相机参数。配置中的 `layerImageDir` 相对仓库根目录。
- `src/` 为主脚本和批量接口 `runUDCImagingV2`，`src/private/` 为透射膜堆求解。
- `Output/<项目>/<模式>/` 保存逐波长和通道透过率、水平 MTF、含二维 PSF/MTF 和线性图像的 MAT 文件，以及无遮挡、同曝光、亮度匹配三张预览。
- `dist/` 为交付预留。运行时需要同仓库的 `topic_coe_reflection` 和 `topic_udc_imaging`，分别提供共用光学算法、D65 数据及线性图像卷积。

透射掩膜只接受纯黑和纯白，分别对应归一化的 0 和 1。普通 8 位图片的纯白编码为 255。图像中出现中间灰阶或彩色像素时，读取阶段直接报错。黑区复振幅透射系数为 0。白区按其光学配置计算，`invert` 固定为 `false`。

- `method = "transmittance"` 用 `transmittance` 指定白区的强度透过率 (T)，取 0 到 1 的标量或逐波长数组。复振幅系数取 `sqrt(T)`，材料相位为零。直接透过率描述整个局部区域，白区不可与另一透射层重叠。
- `method = "nk"` 用 `n`、`k`（标量或逐波长）和 `thickness`（m）设置白区膜层。结合同一配置中的 `incidentN` 和 `exitN`，按从下到上的层序递推法向复透射，计入界面、多次反射和膜厚相位。所有掩膜都为黑的像素保持遮光。

薄膜材料也可沿用 COE 的目录和引用方式，将材料光谱放在项目的 `Materials/` 下，在图层里设 `"nkFile": "Materials/BM_nk.csv"`。路径相对于图层图片目录（即配置中的 `layerImageDir`），膜厚仍通过该图层的 `thickness` 设置。CSV 可以直接写三列数值，依次是波长（nm）、`n`、`k`。波长须严格递增并覆盖整个计算波段。程序进行线性插值，超出表格范围时报错。非空的 `nkFile` 优先于内联值。 [case01 示例材料表](../Input/case01/Materials/BM_nk.csv)

两种方式都可叠加 `heightMap` 指定的单程相位。共用 `displayoptics.readLayerGeometry` 完成采样、整周期裁取或平铺，并生成膜层覆盖区域。

两条流程均以指定均匀频率网格上的傅里叶积分和接收格内中点细分积分计算功率。V2 把接收格定在传感器坐标上，按 `传感器坐标 / 焦距` 转换为近轴角度。通道核按 D65 光谱和相机响应加权，保留透过率及画幅截断损失。单位和 PSF 用于亮度匹配，MTF 从该 PSF 求得。

`case01/transmission_mask.png` 按原 COE 版图生成，阳极区为黑，BM 覆盖区为白，示例在 BM 区验证非零透射。直接模式取 `T=0.10`，有限膜模式取 `n=1.5`、`k=0.2`、厚度 0.5 μm。两种模式均使用 8×8 Pitch（总计 64 Pitch）、64 点/Pitch 和 D65 的 380–780 nm（间隔 10 nm）。相机取 0.20 mm 入瞳直径、3 mm 焦距、1 μm 像元及 1201×1201 接收格。像元中心覆盖 ±0.6 mm，与 Yang 等图 4 的 1.2 mm PSF 显示宽度对齐。RGB 响应采用示例高斯权重。相机和透射材料值均为流程验证参数。

## 阵列与入瞳对比

[MATLAB 算例](../benchmarks/array_pupil/run_comparison.m) 用同一份 `case01` 的 Mask、`nk` 膜层和相机配置运行三组试验。为复现下表，算例固定使用原来的 401×401 接收格。

| 阵列 | 入瞳直径 | G 通道 PSF 水平半高全宽 | G 通道水平 MTF50 |
|---|---:|---:|---:|
| 4×4 | 0.20 mm | 8.69 μm | 44.8 cycles/mm |
| 8×8 | 0.20 mm | 8.69 μm | 44.8 cycles/mm |
| 8×8 | 0.40 mm | 4.32 μm | 92.8 cycles/mm |

从算例脚本生成 MAT 文件后，用 Python 运行 [绘图脚本](../benchmarks/array_pupil/plot_comparison.py)（需要 NumPy、SciPy、Matplotlib），在 `Output/case01/array_pupil_comparison/` 生成 `comparison.png` 和 `summary.csv`。图中的场景按通道匹配亮度。相对信号用窗口捕获比例乘以采样入瞳面积，按相同入射辐照度和曝光估算。

## 成像质量评价

主脚本按相同入瞳、焦距、光谱和传感器网格计算含 Mask 与无遮挡参考。每个模式的输出包括以下资料。

| 文件 | 内容 |
|---|---|
| `udc_channel_metrics.csv` | 各通道总透过率、窗口捕获透过率、无遮挡窗口捕获率、相对无遮挡的捕获能量及窗口覆盖透射能量的比例。 |
| `udc_mtf_directional_summary.csv`、`udc_mtf_directional.csv` | 0°、90°、45°、135° 的 MTF50、频段最小 MTF 和含 Mask／无遮挡的采样曲线。频率沿截线的物理长度计，单位 cycles/mm。MTF50 取从零频开始首次下降到 0.5 的位置。 |
| `udc_psf_tail.csv` | 各半径之外、传感器计算窗以内的绝对尾部透过率与能量占比。占比以各自窗口内捕获能量为分母。 |
| `udc_diffraction_diagnostics.png` | 有效振幅入瞳、RGB 衍射强度、G 通道对数 PSF、G 通道二维 MTF 四联图。 |
| `udc_hdr_point.png` | 黑背景上的无遮挡／含 Mask HDR 点光源，使用共同曝光、空间范围和对数显示映射。 |
| `udc_imaging_comparison.png` | 无遮挡／含 Mask 的 G 通道 PSF 全范围、中心放大和四方向最低 MTF 曲线。PSF 使用共同色标。 |
| `udc_results.mat` | `optics.quality` 保留参考 PSF、参考 MTF、质量表格和点光源的未裁剪线性数据。 |

`cfg.quality.tailRadiiUm`、`mtfBandCyclesPerMm` 和 `hdrPointSignal` 默认分别为 `[10,20,50]`、`[10,100]` 和 `1000`。半径或频段超出有效计算范围时，对应统计记为 NaN。未找到 MTF50 交点时也记为 NaN。跨 Mask 比较应统一半径、频段和计算窗口。窗口捕获比例有助于判断当前窗口是否覆盖足够长的衍射尾部。

点光源是输入线性强度为 `hdrPointSignal` 的理想居中脉冲，响应等于该强度乘以保留能量的衍射核。点光源预览只用于显示动态范围，线性输出未加入传感器饱和、散射或 ISP。相关论文与图示的对应关系见 [文献与评价方法](UDC_diffraction_references.md)。

## 自动 HTML 报告

主脚本完成所有选定模式后，在 `Output/<项目>/udc_report.html` 生成一份汇总报告。每种模式先并排展示无遮挡和含 Mask 的同曝光成像，并对比 PSF 全范围、中心放大和 MTF。报告同时包含三类评价指标、衍射四联图、HDR 点光源和亮度匹配成像。图片嵌入 HTML，单独打开该文件即可离线查看。

每个模式都会计算同一入瞳、焦距、光谱和像元网格下的无遮挡能量核。项目主流程用同一张场景图生成 `udc_clear.png` 和含 Mask 的同曝光图。直接调用 `runUDCImagingV2` 时，提供场景图即可同时得到 `images.clear`。

如果已有完整的 CSV、MAT 和 PNG，可将本模块的 `src` 加入 MATLAB 路径，再调用 `writeUDCReport(项目输出目录, ["specifiedTransmittance","nk"])` 单独重新导出报告。历史 `direct` 目录也可用于生成报告。

## 4 mm 入瞳对比

将本模块的 `src` 加入 MATLAB 路径后，可分别运行 `runUDCProject("case1")` 至 `runUDCProject("case4")`。已有 Case1–3 的同参数结果时，只需运行 `runUDCProject("case4")`。四组共用 236.2 μm Pitch、4 mm 入瞳、3 mm 焦距、1 μm 像元、201×201 接收格、D65 380–780 nm 的 41 个波长和 USAF500 场景。原图每 Pitch 1000 像素，`cfg.samplesPerPeriod=500` 按像素中心最近邻降采样，再平铺到 17×17 Pitch。HDR 点光源预览通过 `cfg.quality.hdrDisplayHalfWidthUm=50` 单独放大中央区域，指标和线性数据仍以完整计算窗口为准。

- Case1 与 Case3 的黑区为 2 μm BM，白区没有 BM。Case2 每张图的黑区对应该膜层，12 nm Cath 在下，2 μm BM 在上。均设 `stack.backgroundMode="medium"`，无膜区按两侧介质界面计算透射。
- Case4 使用 Case1 的 BM 图样及材料，沿用 Case2 的阴极材料、厚度和入射／出射介质。`coverage="full"` 在 BM 下方的每个采样点设置阴极，覆盖全部物面及 BM 开口。`Input/case4/cath.png` 为保留的原图。
- `materialOn="black"` 选择原始二值图的黑区覆盖材料。默认图层仍用白区标记材料覆盖。材料表以实际波长节点线性插值，可读取逗号或空白分隔的三列数据。
- 大入瞳算例设置 `cfg.receiverMethod="fftBins"`，分块预计算入瞳各区域落在接收窗口内的傅里叶基底，再依据 FFT 原生格与传感器像元的重叠分配能量。`cfg.receiverFftPaddingFactor` 控制零填充，默认的中点求积使用 `receiverMethod="midpoint"`。4 mm 入瞳结果属于当前标量、近轴、共面模型的方案对比。

四组结果齐备后，可用 Python 运行[对比绘图脚本](../benchmarks/case_comparison/plot_comparison.py)（需要 h5py、NumPy、Matplotlib、Pillow）。脚本验证共同的无遮挡参考和 Case4 的实际覆盖，并在 `Output/case1_case2_case3_case4_4mm_500pitch/` 生成 `comparison.png` 和 `summary.csv`。原三组汇总仍可在 `Output/case1_case2_case3_4mm_500pitch/` 查阅。

## 半径 59 μm 与 54 μm 的开口对比

运行 [Case5–8 输入生成脚本](../benchmarks/case_comparison/prepare_aperture_cases.py) 可复现输入。脚本用 Case1 二值图的中心圆和周期边界圆拟合圆心，在原始像素中心采样 54 μm 圆。Case5 的 59 μm BM 沿用 Case1 原图，54 μm 圆为 Cath 开口。Case6 的 54 μm 圆为 BM 开口，膜堆只有 BM。

依次运行 `runUDCProject("case5")`、`runUDCProject("case6")`、`runUDCProject("case7")` 和 `runUDCProject("case8")`。Case7 引用已计算的 Case4，配置和膜层输入相同。Case8 引用 Case7 并乘以理想中性衰减，场振幅乘 √0.5，绝对能量和同曝光线性图像乘 0.5。归一化 PSF、MTF、亮度匹配图像和无遮挡参考沿用基准。

给[对比绘图脚本](../benchmarks/case_comparison/plot_comparison.py)传入 `--apertures`，生成 `Output/case1_case5_case6_case7_case8_4mm_500pitch/` 下的对比图和包含无遮挡行的汇总 CSV。脚本核对 54 μm 掩膜降采样后的覆盖、连续阴极覆盖和中性衰减比例。

## 交互式结果对比

已生成的 UDC V2 结果可在同级的[透射衍射对比工作台](../../topic_transmission_diffraction_viewer/docs/README.md)中查看、打包和比较。工作台会自动发现本 topic 的 `Output/` 中的结果。
