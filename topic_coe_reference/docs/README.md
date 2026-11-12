# 原作者 COE 参考实现（暂停开发）

本方向目前暂停开发，保留已有代码与资料供功能设计参考。后续主要开发方向为 `topic_coe_reflection` 和 `topic_udc_imaging`。

本目录为同事编写的像素显示仿真代码，主要借鉴其功能设计和使用流程。后续 COE/UDC 开发的计算算法按本项目的物理模型与验证要求独立评估和实现。

此目录保留原 `WangAnPing_code/衍射` 的两版程序、材料文本和原始压缩包；像素展开 GUI 已独立移到 `topic_pixel_unfolding`。

两版文件改为合法且不同的 MATLAB 函数名：

```matlab
addpath('topic_coe_reference/src');
CoE_Diffraction_v8_spectrum_1;
% 或 CoE_Diffraction_v8_spectrum_2;
```

各文件的计算辅助函数均包含在文件内部；输入从本方向 `Input/` 读取，输出写入本方向 `Output/`，与工作目录无关。物理公式保留原作者实现；维护版反射计算请用 `topic_coe_reflection`。

五张掩模 `BM.png`、`PDL.png`、`RCF.png`、`GCF.png`、`BCF.png` 已从原始压缩包提取到 `Input/`，材料文本也集中到该目录。原程序还要求 `RCF.txt`、`GCF.txt` 和 `BCF.txt`。现有红色材料数据名为 `RCF02.txt`，压缩包中同样没有 `RCF.txt`；不能无依据将两者视为同一材料。运行参考实现前须确认并补齐 `RCF.txt`。蓝色材料读取路径已与现有 `bcf.txt` 的大小写一致。

`dist/衍射.7z` 原样保留。旧目录中剩余的原始掩模图片保存在 `Input/original_source/衍射/`，避免覆盖已整理的输入。此目录不随 `setupDiffraction` 自动加入路径，也不作为本次光学回归验证对象。原图像处理/弱信号显示流程需要 Image Processing Toolbox。

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
