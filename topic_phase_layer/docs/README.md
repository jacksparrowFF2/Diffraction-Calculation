# 相位层与 DOE 优化（暂停开发）

本模块是此前为专利工作开发的优化程序，目前暂停开发，保留已有功能供查阅和按需复用。后续主要开发方向为 `topic_coe_reflection` 和 `topic_udc_imaging`。

```matlab
addpath('topic_phase_layer/src');
CustomApertureWithPhase_fraunhofer_fft;
```

入口默认尝试读取本方向 `Output/` 中的 DOE 结果；没有结果时使用脚本内参数。`useBestDoeDesign=false` 可直接用手动设计。

```matlab
OptimizeMultiRadiusPhaseIslands; % 默认 full 全量搜索，耗时可能较长
design = loadBestMultiRadiusPhaseDesign();
```

首次调试建议将优化脚本开头的 `doePreset` 改为 `"quick"`。`src/private/` 自带原孔径/传播/卷积依赖，`Input/` 带默认图，`Output/` 保存 DOE CSV/MAT 与可视化；整目录可独立复制。物理模型和原优化目标保持既有行为。

历史 DOE 结果已迁移到本方向 `Output/`。读取其它位置的结果时，可显式指定路径：

```matlab
design = loadBestMultiRadiusPhaseDesign('旧结果文件或目录的完整路径');
```

图像处理功能需要 Image Processing Toolbox。完整 DOE 不在自动测试中执行；自动检查覆盖功能目录独立运行和结果读取路径。

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
