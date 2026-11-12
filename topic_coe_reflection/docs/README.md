# COE 反射衍射

本模块是后续主要开发方向之一，用于评估产品 COE 反射设计方案的可行性，并在一致的计算条件与评价指标下对比不同方案的效果。

本目录可独立复制，只需添加本方向 `src/` 路径，`src/private/` 随目录一起保留。

```matlab
addpath('topic_coe_reflection/src'); % 在仓库根目录执行；独立复制后改成实际路径
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
