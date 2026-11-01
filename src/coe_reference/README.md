# 原作者 COE 参考实现

此目录保留原 `WangAnPing_code/衍射` 的两版程序、材料文本和原始压缩包；像素展开 GUI 已独立移到 `pixel_unfolding`。

两版文件改为合法且不同的 MATLAB 函数名：

```matlab
addpath('src/coe_reference');
CoE_Diffraction_v8_spectrum_1;
% 或 CoE_Diffraction_v8_spectrum_2;
```

各文件的计算辅助函数均包含在文件内部；输入从本目录 `input/` 读取，输出写入本目录 `output/`，与工作目录无关。物理公式保留原作者实现；维护版反射计算请用 `coe_reflection`。

五张掩模 `BM.png`、`PDL.png`、`RCF.png`、`GCF.png`、`BCF.png` 已从原始压缩包提取到 `input/`，材料文本也集中到该目录。原程序还要求 `RCF.txt`、`GCF.txt` 和 `BCF.txt`。现有红色材料数据名为 `RCF02.txt`，压缩包中同样没有 `RCF.txt`；不能无依据将两者视为同一材料。运行参考实现前须确认并补齐 `RCF.txt`。蓝色材料读取路径已与现有 `bcf.txt` 的大小写一致。

`archive/衍射.7z` 原样保留。此目录不随 `setupDiffraction` 自动加入路径，也不作为本次光学回归验证对象。原图像处理/弱信号显示流程需要 Image Processing Toolbox。
