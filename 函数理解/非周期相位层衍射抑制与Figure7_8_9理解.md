# 非周期相位层衍射抑制与 Figure 7/8/9 理解

本文档用于记录 `CustomApertureWithPhase_fraunhofer_fft.m` 中非周期相位去相关层仿真的判读思路，便于后续继续讨论和修改代码。

## 1. 抑制衍射的目标

对于 COE OLED 息屏彩虹纹、彩色衍射条纹、点光源反射星芒等问题，所谓“抑制衍射”不是把衍射能量集中到某一个新的峰，也不是让总能量消失。

结合实际显示/成像效果，更准确的目标应分主次：

```text
优先把有效能量保留在 0 级主斑/镜面主反射附近；
减少可见的非零级衍射环、星芒和规则条纹数量；
降低特定衍射级次峰值；
降低非零级衍射峰相对 0 级主斑的强度比例；
降低 RGB 分离造成的彩色亮点或彩纹；
必要时才把无法消除的尖锐、方向性强的非零级衍射峰展宽并打散为低强度背景；
同时避免过强宽角散射导致 haze、sparkle 或清晰度下降。
```

因此，相位去相关层的理想效果不是单纯“把能量分散开”，而是“0 级集中、非零级削弱、少衍射环、低 haze”。展宽/去相关只是对非零级强峰的补救手段，不应以牺牲 0 级能量和清晰度为代价。

换句话说，如果两个方案都能压低彩色衍射环，那么衍射环更少、0 级主斑能量更集中、宽角背景更低的方案更接近实际最优；不是红色差分区域越多、能量越分散越好。

如果相位结构本身具有规则周期，可能会产生新的衍射级次；这不是本方案希望看到的结果。因此当前代码采用泊松盘非周期圆形高折微岛，而不是规则微透镜阵列、规则光栅或规则圆柱阵列。

## 2. 为什么 before/after 肉眼不一定明显

当前远场图中，基础孔径 `Mask` 仍然是主要的周期性振幅结构。非周期相位层通过如下复振幅项叠加：

```matlab
ComplexMask = BaseMask .* exp(1j * phi_double);
```

它会调制相干峰并改变能量分布，但不会完全移除原始孔径周期带来的远场骨架。实际判读时不能只看非零级峰是否被削弱，还要看 0 级主斑是否被明显削弱、变宽或把能量推到宽角背景。

此外，普通 RGB 图或强度图常使用归一化、对数压缩和伽马增强。这些显示方法有利于看清图案，但也可能把峰值降低的差异压缩掉。例如绿色通道当前输出为：

```text
550 nm: 无相位 0.0001302, 有相位 0.0001017, 峰值比例 78.07 %
```

这说明峰值确实降低，但在普通图像显示中未必很醒目。因此需要 Figure 7、8、9 这类诊断图。

## 3. Figure 7：Green channel before/after 统一 dB 色标

对应导出文件：

```text
ImageForShow/green_before_after_unified_db.png
```

这张图左右分别为：

```text
左图：无相位层绿色通道远场强度；
右图：有非周期相位层绿色通道远场强度。
```

两张图使用同一个参考最大值和同一个 `[-60, 0] dB` 色标，因此可以直接比较亮度高低。

判读重点：

```text
原来的亮峰是否变暗；
尖锐峰是否变钝或展宽；
0 级主斑是否仍然集中；
峰周围是否出现过宽的低强度背景；
是否产生新的强规则峰。
```

这张图回答的问题是：

```text
相位层有没有压低非零级高强度尖峰，同时尽量保持 0 级主斑集中，并避免形成明显宽角背景。
```

## 4. Figure 8：Green channel 差分图

对应导出文件：

```text
ImageForShow/green_phase_difference.png
```

差分定义为：

```matlab
GreenDiff = (GreenIntensity_WithPhase - GreenIntensity_NoPhase) / greenReference;
```

颜色含义：

```text
蓝色：加相位后该角度能量降低；
白色：变化较小；
红色：加相位后该角度能量增强。
```

如果在原非零级衍射峰附近看到蓝色，而峰周围或其他角度看到红色，说明相位层把原本集中在非零级衍射峰中的能量转移到了更宽的角度范围。

这张图回答的问题是：

```text
相位层把能量从哪里削弱，又把能量重新分布到了哪里。
```

需要注意：红色区域不是自动等于好事。若红色主要出现在 0 级主斑附近，且表示能量回到主斑或主反射方向，通常是有利的；若红色形成新的尖锐强峰，或在宽角区域形成明显背景，则可能带来新的衍射、haze 或清晰度问题。对实际方案而言，应优先接受“非零级变蓝、0 级保持或增强、宽角红色少”的结果。

## 5. Figure 9：中心水平截面强度曲线

对应导出文件：

```text
ImageForShow/green_center_horizontal_section_linear_db.png
```

这张图取绿色通道远场图的中心水平线，比较无相位层和有相位层的强度曲线。

曲线含义：

```text
黑线：无相位层；
红线：有非周期相位层。
```

上半图为线性归一化强度，适合观察主峰和较强峰的高度变化。下半图为 dB 强度，适合观察弱旁瓣和低强度背景变化。

判读重点：

```text
红线在某些峰位低于黑线：该衍射峰被削弱；
红线在 0 级主斑附近高于或接近黑线：主斑能量保留较好；
红线在峰周边或远离峰的位置高于黑线：能量被展宽或重新分布，需要结合 haze 风险判断；
红线如果在新位置形成尖锐高峰：可能引入新的衍射问题，需要避免。
```

这张图回答的问题是：

```text
沿一个固定方向，相位层具体削弱了哪些峰，又增强了哪些背景或旁瓣。
```

## 6. 当前设计是否最优

当前泊松盘圆形高折微岛设计不是最优设计，而是第一版 baseline：

```text
可量产友好；
参数少；
可复现；
物理含义清楚；
便于作为专利概念验证和 DOE 起点。
```

它的局限包括：

```text
只有单一半径；
只有单一高度；
只有二值相位；
没有显式优化空间频谱；
没有针对像素/BM/touch 主频避让；
没有针对 RGB 三波长共同优化；
没有针对目标衍射峰做反向优化。
```

后续更优方向可以包括：

```text
多半径圆形微岛；
多高度灰阶相位岛；
随机旋角椭圆微岛；
短条形或短弧形相位结构；
带通/蓝噪声频谱约束高度图；
以“非零级目标衍射峰最小 + 0 级主斑能量保持 + haze 风险受限”为目标函数的优化设计。
```

## 7. DOE 最优设计如何导入可视化程序

`OptimizeMultiRadiusPhaseIslands.m` 会把 DOE 结果保存到：

```text
ImageForShow/multi_radius_phase_doe_results.csv
ImageForShow/multi_radius_phase_doe_results.mat
```

CSV 中每一行代表一个候选相位层设计。当前脚本最后执行：

```matlab
results = sortrows(results, 'Objective', 'ascend');
```

因此，CSV 文件中 `Objective` 最小的行就是当前目标函数下的最优设计；由于已经按 `Objective` 升序排序，通常第一行数据就是当前 DOE 的最优候选。

需要注意：这里的“最优”是当前 DOE 搜索范围、当前目标函数权重、当前 haze 惩罚设定下的最优，不等同于全局物理最优。

CSV 字段与可视化脚本参数的对应关系如下：

```text
RadiusSet_um        -> phaseIslandRadiusList
TargetFillFactor    -> phaseFillFactor
MinDistanceFactor   -> phaseMinDistanceFactor
Seed                -> phaseSeed
```

旧版 quick DOE 的候选示例为：

```text
RadiusSet_um = [1 1.5]
TargetFillFactor = 0.40
MinDistanceFactor = 1.00
Seed = 1
ActualFillFactor ≈ 0.3818
TopKPeakRatio ≈ 0.6463
HazeRisk ≈ 0.1489
Objective ≈ 0.7210
```

新版 DOE 会额外输出 `NonZeroTopKPeakRatio`、`ZeroOrderRetention` 和 `ZeroOrderLossPenalty`。当前 `CustomApertureWithPhase_fraunhofer_fft.m` 已可通过 `loadBestMultiRadiusPhaseDesign.m` 自动导入最优候选；如需手动设置，则对应关系为：

```matlab
phaseIslandRadiusList = [1.0, 1.5] * 1e-6;
phaseFillFactor = 0.40;
phaseMinDistanceFactor = 1.00;
phaseMinDistance = phaseMinDistanceFactor * 2 * max(phaseIslandRadiusList);
phaseSeed = 1;
```

然后重新运行可视化脚本，即可生成该候选设计对应的相位 map、3D 相位分布、Figure 7/8/9 以及峰值强度对比。

还需要特别关注 `HazeRisk` 和 0 级主斑保持情况。如果某个候选设计 `Objective` 排名靠前，但 `HazeRisk` 明显超过 `hazeRiskLimit`，或 0 级主斑能量明显下降，说明它可能属于“非零级削峰强，但把能量打散过多”的方案。此时不能只看 `Objective`，还应结合 0 级能量、haze 风险和 Figure 7/8/9 判断是否适合作为实际设计。

## 8. 当前结论

当前仿真中，非周期相位层的作用可以概括为：

```text
不是消灭衍射图案；
不是把能量集中到新峰；
也不是单纯把能量分散得越开越好；
而是尽量保持 0 级主斑能量集中，同时降低非零级衍射环/条纹/星芒的可见度，并限制宽角散射背景。
```

后续看图时，应优先结合：

```text
Figure 7：整体 before/after dB 对比；
Figure 8：能量变化位置；
Figure 9：中心截面定量曲线；
0 级主斑能量保持率；
非零级峰值强度比例和非零级峰值/背景比；
haze 或宽角散射风险。
```
