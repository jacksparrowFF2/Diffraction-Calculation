# UDC 透射衍射成像 V2

在 MATLAB 中运行 [UDCImagingV2_main.m](../src/UDCImagingV2_main.m)。默认读取 [case01 配置](../Input/case01/config.json)，以同一张透射掩膜和 V1 的 `USAF500.png` 场景运行 `direct` 和 `nk` 两种配置。主脚本的 `transmissionModes` 可只选择其中一种。配套运行 [COEReflection_main.m](../../topic_coe_reflection/src/COEReflection_main.m) 可得到同一版图的反射衍射。

- `Input/<项目>/` 存放透射掩膜和相机参数。配置中的 `layerImageDir` 相对仓库根目录。
- `src/` 为主脚本和批量接口 `runUDCImagingV2`，`src/private/` 为透射膜堆求解。
- `Output/<项目>/<模式>/` 保存逐波长和通道透过率、水平 MTF、含二维 PSF/MTF 和线性图像的 MAT 文件，以及无遮挡、同曝光、亮度匹配三张预览。
- `dist/` 为交付预留。运行时需要同仓库的 `topic_coe_reflection` 和 `topic_udc_imaging`，分别提供共用光学算法、D65 数据及线性图像卷积。

透射掩膜只接受纯黑和纯白，分别对应归一化的 0 和 1。普通 8 位图片的纯白编码为 255。图像中出现中间灰阶或彩色像素时，读取阶段直接报错。黑区复振幅透射系数为 0。白区按其光学配置计算，`invert` 固定为 `false`。

- `method = "transmittance"` 用 `transmittance` 指定白区的强度透过率 (T)，取 0 到 1 的标量或逐波长数组。复振幅系数取 `sqrt(T)`，材料相位为零。直接透过率描述整个局部区域，白区不可与另一透射层重叠。
- `method = "nk"` 用 `n`、`k`（标量或逐波长）和 `thickness`（m）设置白区膜层。结合同一配置中的 `incidentN` 和 `exitN`，按从下到上的层序递推法向复透射，计入界面、多次反射和膜厚相位。所有掩膜都为黑的像素保持遮光。

两种方式都可叠加 `heightMap` 指定的单程相位。共用 `displayoptics.readLayerGeometry` 完成采样、整周期裁取或平铺，并生成膜层覆盖区域。

两条流程均以指定均匀频率网格上的傅里叶积分和接收格内中点细分积分计算功率。V2 把接收格定在传感器坐标上，按 `传感器坐标 / 焦距` 转换为近轴角度。通道核按 D65 光谱和相机响应加权，保留透过率及画幅截断损失。单位和 PSF 用于亮度匹配，MTF 从该 PSF 求得。

`case01/transmission_mask.png` 按原 COE 版图生成，阳极区为黑，BM 覆盖区为白，示例在 BM 区验证非零透射。直接模式取 `T=0.10`，有限膜模式取 `n=1.5`、`k=0.2`、厚度 0.5 μm。两种模式均使用 4×4 Pitch、64 点/Pitch 和 D65 的 380–780 nm（间隔 10 nm）。相机取 0.20 mm 入瞳直径、3 mm 焦距、1 μm 像元及 401×401 接收格，RGB 响应采用示例高斯权重。相机和透射材料值均为流程验证参数。
