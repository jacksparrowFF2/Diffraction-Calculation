# UDC 成像

本目录可独立复制；保留 `private/` 和 `input/`，只添加本目录路径。

```matlab
addpath('src/udc_imaging');
cfg = udcImagingDefaults();
[optics, images] = runUDCImaging(cfg); % 本目录 input/USAF500.png
image(images.fixedExposure); axis image;

optics = computeUDCOptics(cfg); % 只计算 PSF/MTF
hdr = zeros(501,501,3);
hdr(251,251,:) = 1000;
images = formUDCImage(hdr, optics, "linear");
```

普通图片使用 `runUDCImaging(cfg, imread('scene.png'), "srgb")`。函数返回线性数组，不自动保存；`fixedExposure` 保留损耗，`brightnessMatched` 仅比较单位和 PSF 模糊。亮度补偿结果不包含传感器噪声或恢复算法。

`private/` 自带计算核心、默认参数、孔径生成和 FFT 卷积，无需其他功能目录。采用法向标量、屏幕与入瞳共面、远景对焦模型；默认参数是示例值。完整说明见仓库 `docs/COE反射与UDC成像分支.md`。
