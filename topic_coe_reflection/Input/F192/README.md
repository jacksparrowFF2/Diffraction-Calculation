# F192：实际版图输入

- 主脚本设置 `projectName="F192"`，保持 `useProjectConfig=true`。
- 六张 PNG 均为 1000×1000，覆盖每轴 4 Pitch；每 Pitch 250 像素。
- config 中 Pitch 为 59.06 μm，当前目标 `[32,32]`、`repeatToFit=true`：横纵各重复8次，计算网格8000×8000；不缩放或改写原图。
- anode、CF_R/B/G 白色表示材料；BM、IJP 黑色表示材料，配置均为 `invert=true`。全黑 IJP 反相后整面覆盖。
- CF 从下到上为 R→B→G，各自占独立层级；交叠处保留各层完整膜厚。自定义顺序、光学假设及输入规则见 [统一使用说明](../../docs/README.md)。
- n/k 和厚度目前仍为配置中的占位参数；Materials 中的数据仅在填写相应 `nkFile` 后启用。

- 当前保留250像素/Pitch；32Pitch用于D65固定ROI方案比较，相对原64Pitch通过宽带容差，不能宣称逐波长或峰形已收敛。64Pitch可作为复核设置，实际材料改变后须重跑收敛基准。
