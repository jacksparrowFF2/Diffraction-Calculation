# COE 反射衍射

本目录可独立复制，只需添加本目录路径，`private/` 随目录一起保留。

```matlab
addpath('src/coe_reflection'); % 在仓库根目录执行；独立复制后改成实际路径
cfg = coeReflectionDefaults();
cfg.featureIntensity = 0.60;
cfg.backgroundIntensity = 0.01;
r = runCOEReflection(cfg);
disp(r.metrics);
imagesc(rad2deg(r.coordinate), rad2deg(r.coordinate), r.kernel(:,:,2));
axis image; colorbar;
```

`heightMap` 以米为单位，反射使用双程相位；支持自定义 `mask` 和复数 `complexCoefficient`。`baselineKernel/baselineMetrics` 为清除额外高度图后的基准。函数返回物理数组，不自动写文件；需要时可用 `save` 保存。

本目录自带默认参数、FFT/能量积分核心和几何生成代码，不依赖 UDC、旧孔径或根目录辅助函数。法向、标量、近轴远场简化；默认光谱与反射率为示例值。完整指标解释见仓库 `docs/COE反射与UDC成像分支.md`。
