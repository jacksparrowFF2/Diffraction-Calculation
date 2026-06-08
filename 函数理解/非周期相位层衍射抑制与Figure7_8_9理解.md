# 非周期相位层衍射抑制与 Figure 7/8/9 理解

本文档用于记录 `CustomApertureWithPhase_fraunhofer_fft.m` 中非周期相位去相关层仿真的判读思路，便于后续继续讨论和修改代码。

## 1. 抑制衍射的目标

对于 COE OLED 息屏彩虹纹、彩色衍射条纹、点光源反射星芒等问题，所谓“抑制衍射”不是把衍射能量集中到某一个新的峰，也不是让总能量消失。

更准确的目标是：

```text
降低特定衍射级次峰值；
降低峰值/背景比；
降低 RGB 分离造成的彩色亮点或彩纹；
把尖锐、方向性强的衍射峰展宽并打散为低强度背景；
同时避免过强宽角散射导致 haze、sparkle 或清晰度下降。
```

因此，相位去相关层的理想效果是“削峰、展宽、去相关、低 haze”，而不是把能量导流到另一个更亮、更规则的新衍射峰。

如果相位结构本身具有规则周期，可能会产生新的衍射级次；这不是本方案希望看到的结果。因此当前代码采用泊松盘非周期圆形高折微岛，而不是规则微透镜阵列、规则光栅或规则圆柱阵列。

## 2. 为什么 before/after 肉眼不一定明显

当前远场图中，基础孔径 `Mask` 仍然是主要的周期性振幅结构。非周期相位层通过如下复振幅项叠加：

```matlab
ComplexMask = BaseMask .* exp(1j * phi_double);
```

它会削弱相干峰、展宽部分能量，但不会完全移除原始孔径周期带来的远场骨架。

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
峰周围是否出现更宽的低强度背景；
是否产生新的强规则峰。
```

这张图回答的问题是：

```text
相位层有没有把远场中的高强度尖峰削弱，并让能量分布变得更分散。
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

如果在原衍射峰附近看到蓝色，而峰周围或其他角度看到红色，说明相位层把原本集中在衍射峰中的能量转移到了更宽的角度范围。

这张图回答的问题是：

```text
相位层把能量从哪里削弱，又把能量重新分布到了哪里。
```

需要注意：红色区域并不一定是坏事。只要红色不是形成新的尖锐强峰，而是低强度、宽角度背景，就符合“削峰展宽”的目标。

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
红线在峰周边或远离峰的位置高于黑线：能量被展宽或重新分布；
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
以“目标衍射峰最小 + haze 风险受限”为目标函数的优化设计。
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

例如 quick DOE 当前最优候选为：

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

则需要在 `CustomApertureWithPhase_fraunhofer_fft.m` 中设置：

```matlab
phaseIslandRadiusList = [1.0, 1.5] * 1e-6;
phaseFillFactor = 0.40;
phaseMinDistanceFactor = 1.00;
phaseMinDistance = phaseMinDistanceFactor * 2 * max(phaseIslandRadiusList);
phaseSeed = 1;
```

然后重新运行可视化脚本，即可生成该候选设计对应的相位 map、3D 相位分布、Figure 7/8/9 以及峰值强度对比。

还需要特别关注 `HazeRisk`。如果某个候选设计 `Objective` 排名靠前，但 `HazeRisk` 明显超过 `hazeRiskLimit`，说明它可能属于“削峰强，但 haze 风险偏高”的方案。此时不能只看 `Objective`，还应结合 haze 风险和 Figure 7/8/9 判断是否适合作为实际设计。

## 8. 当前结论

当前仿真中，非周期相位层的作用可以概括为：

```text
不是消灭衍射图案；
不是把能量集中到新峰；
而是降低原有强峰、展宽部分能量、降低方向性衍射峰的可见度。
```

后续看图时，应优先结合：

```text
Figure 7：整体 before/after dB 对比；
Figure 8：能量变化位置；
Figure 9：中心截面定量曲线；
峰值强度比例和峰值/背景比等数值指标。
```
