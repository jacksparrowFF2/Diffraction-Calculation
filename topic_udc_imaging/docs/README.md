# UDC 成像

本模块是后续主要开发方向之一，用于评估产品 UDC 成像设计方案的可行性，并在一致的场景、曝光条件与评价指标下对比不同方案的效果。

本目录可独立复制；保留 `src/private/` 和 `Input/`，只添加本方向 `src/` 路径。

```matlab
addpath('topic_udc_imaging/src');
cfg = udcImagingDefaults();
[optics, images] = runUDCImaging(cfg); % 本方向 Input/USAF500.png
image(images.fixedExposure); axis image;

optics = computeUDCOptics(cfg); % 只计算 PSF/MTF
hdr = zeros(501,501,3);
hdr(251,251,:) = 1000;
images = formUDCImage(hdr, optics, "linear");
```

普通图片使用 `runUDCImaging(cfg, imread('scene.png'), "srgb")`。函数返回线性数组，不自动保存；`fixedExposure` 保留损耗，`brightnessMatched` 仅比较单位和 PSF 模糊。亮度补偿结果不包含传感器噪声或恢复算法。

`src/private/` 自带计算核心、默认参数、孔径生成和 FFT 卷积，无需其他功能目录。采用法向标量、屏幕与入瞳共面、远景对焦模型；默认参数是示例值。完整说明见仓库 `docs/COE反射与UDC成像分支.md`。

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
