# 相位层与 DOE 优化

```matlab
addpath('src/phase_layer');
CustomApertureWithPhase_fraunhofer_fft;
```

入口默认尝试读取本目录 `output/` 中的 DOE 结果；没有结果时使用脚本内参数。`useBestDoeDesign=false` 可直接用手动设计。

```matlab
OptimizeMultiRadiusPhaseIslands; % 默认 full 全量搜索，耗时可能较长
design = loadBestMultiRadiusPhaseDesign();
```

首次调试建议将优化脚本开头的 `doePreset` 改为 `"quick"`。`private/` 自带原孔径/传播/卷积依赖，`input/` 带默认图，`output/` 保存 DOE CSV/MAT 与可视化；整目录可独立复制。物理模型和原优化目标保持既有行为。

旧结果如在仓库 `ImageForShow/`，请复制到本目录 `output/`，或者：

```matlab
design = loadBestMultiRadiusPhaseDesign('旧结果文件或目录的完整路径');
```

图像处理功能需要 Image Processing Toolbox。完整 DOE 不在自动测试中执行；自动检查覆盖功能目录独立运行和结果读取路径。
