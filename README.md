# Diffraction Calculation

MATLAB 周期孔径衍射、非周期相位层优化、COE 反射衍射和 UDC 成像。

## 开始运行

在 MATLAB 中切换到本仓库根目录：

```matlab
addpath('src');
setupDiffraction;
demoCOEUDC;
```

完整演示结果写入 `ImageForShow/COE_UDC`。只运行某个功能时，可以只添加该功能目录：

```matlab
addpath('src/coe_reflection');
cfg = coeReflectionDefaults();
reflection = runCOEReflection(cfg);

addpath('src/udc_imaging');
cfg = udcImagingDefaults();
[udc, images] = runUDCImaging(cfg);
```

## 按功能组织的源代码

| 目录 | 功能与入口 |
|---|---|
| [aperture_diffraction](src/aperture_diffraction/README.md) | 原始自定义孔径衍射；`CustomAperture_fraunhofer_fft`、四孔预览和显示映射 |
| [phase_layer](src/phase_layer/README.md) | 相位微岛、相位层前后对照、DOE 优化和最优设计读取 |
| [coe_reflection](src/coe_reflection/README.md) | COE 双程反射、多波长衍射、无额外相位层基准及能量统计 |
| [udc_imaging](src/udc_imaging/README.md) | 单程透射、相机 PSF/MTF、RGB/HDR 成像 |
| [pixel_unfolding](src/pixel_unfolding/README.md) | 像素展开与批量图像仿真 GUI |
| [coe_reference](src/coe_reference/README.md) | 原作者 COE 多波长参考实现及材料数据，单独保留 |

每个功能目录保留其所需代码；内部辅助函数放在本目录的 `private/`，MATLAB 自动解析，无需 `addpath` 私有目录。维护版功能目录可单独复制使用；默认测试图已随对应目录提供。`src` 根目录仅放初始化、联合演示和兼容的通用默认参数入口。`reference/` 保留原始教学参考脚本。

基础孔径与相位层脚本的结果写到各自的 `output/`，与工作目录无关；旧 DOE 结果若在 `ImageForShow/`，请复制到 `src/phase_layer/output/` 或给 `loadBestMultiRadiusPhaseDesign` 传入旧文件的完整路径。输出目录不进入 Git。

新增 COE/UDC 模型在 MATLAB R2024a 验证，仅依赖基础 MATLAB。旧图像脚本和 GUI 使用 Image Processing Toolbox；GUI 的 Excel 导出和第三方参考程序另见各自说明。模型假设和参数解释见 [COE/UDC 使用说明](docs/COE反射与UDC成像分支.md)。

## 验证

```matlab
results = runtests('tests');
assertSuccess(results);
```

包含光学解析回归，以及把功能目录复制到临时位置、单独添加路径后的独立运行检查。测试会恢复 MATLAB 路径和界面设置；旧 DOE 完整搜索耗时较长，测试不执行其全量扫描。
