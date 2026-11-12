# Diffraction Calculation

MATLAB 孔径衍射、COE 反射和 UDC 成像仿真。后续开发以 COE/UDC 产品设计方案的可行性评估和不同方案的效果对比为主。

## 方向与定位

| 方向 | 定位 | 说明 |
|---|---|---|
| `topic_aperture_fraunhofer` | 其它方向衍射计算的基础 | [使用说明](topic_aperture_fraunhofer/docs/README.md) |
| `topic_phase_layer`（暂停开发） | 专利工作中的非周期相位层优化程序 | [使用说明](topic_phase_layer/docs/README.md) |
| `topic_coe_reflection` | 主要开发方向：COE 反射方案评估与对比 | [使用说明](topic_coe_reflection/docs/README.md) |
| `topic_udc_imaging` | 主要开发方向：UDC 成像方案评估与对比 | [使用说明](topic_udc_imaging/docs/README.md) |
| `topic_pixel_unfolding` | 同事编写的像素显示仿真 GUI | [使用说明](topic_pixel_unfolding/docs/README.md) |
| `topic_coe_reference`（暂停开发） | 同事编写的像素显示仿真程序，保留供功能设计参考 | [使用说明](topic_coe_reference/docs/README.md) |

`topic_aperture_fraunhofer` 提供基础模型和方法参考。各方向自带运行依赖，不必在运行时加载该基础方向。`topic_coe_reference` 的借鉴重点是功能设计和使用流程，后续算法按本项目的物理模型与验证要求独立评估和实现。

## 统一目录结构

```text
topic_方向名称/
├── Input/           输入图片序列、项目参数文本及其它计算数据
├── Output/          计算产生的图片、表格、MAT 等结果
├── docs/            使用说明、报告
├── dist/            打包发布文件
└── src/             主程序源码
    └── private/     本方向的辅助函数
```

`private` 是命名和层级上的例外：MATLAB 通过这个固定的小写名称识别私有函数，而且它必须位于调用程序所在目录下面。因此保留 `src/private/`，不改成 `PrivateFunction`，不单独添加其路径。只需添加方向的 `src/`，MATLAB 就会自动解析辅助函数。

整个 `topic_…` 文件夹可以独立复制，输入和输出路径由主程序位置定位，与当前工作目录无关。空目录通过 `.gitkeep` 保留；`Output/` 中的新计算结果默认不进入 Git；迁移前已纳入版本管理的历史结果继续保留。

仓库根目录的 `src/` 保留初始化、联合演示和兼容参数入口；`docs/` 保留跨方向的模型说明。`reference/`、`input/` 和 `ImageForShow/` 保留共享教学参考及尚未归属具体方向的历史资料。

## 开始运行

在 MATLAB 中切换到仓库根目录：

```matlab
addpath('src');
setupDiffraction;
demoCOEUDC;
```

联合演示默认写入 `topic_coe_reflection/Output/COE_UDC/`。只运行某个方向时，只添加它的 `src/`：

```matlab
addpath('topic_coe_reflection/src');
cfg = coeReflectionDefaults();
reflection = runCOEReflection(cfg);

addpath('topic_udc_imaging/src');
cfg = udcImagingDefaults();
[udc, images] = runUDCImaging(cfg);
```

独立复制后，将上述路径改为所复制方向的实际 `src/` 路径即可。`setupDiffraction` 不自动加载同事的 COE 参考程序。

相位层入口默认从 `topic_phase_layer/Output/` 读取 DOE 结果；也可给 `loadBestMultiRadiusPhaseDesign` 传入结果文件或目录的完整路径。原 `ImageForShow/` 中能明确归属孔径与相位层方向的结果已迁入各自的 `Output/`。

COE/UDC 计算仅依赖基础 MATLAB；旧图像脚本和 GUI 使用 Image Processing Toolbox。模型假设和参数解释见 [COE/UDC 使用说明](docs/COE反射与UDC成像分支.md)。默认材料系数和光谱均为示例参数，具体产品评估需使用相应设计数据并检查采样收敛。

## 验证

```matlab
results = runtests('tests');
assertSuccess(results);
```

测试包含物理解析回归，以及将整个方向复制到临时位置、只添加其 `src/` 后的独立运行检查。测试恢复路径和界面设置，不执行耗时的完整 DOE 搜索。
