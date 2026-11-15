# Case1_Pitch128：完整128×128 Pitch图片项目

- 来源：Case1六张1024×1024 PNG原样按8×8平铺，未缩放或修改像素值。
- 六张图片均为8192行×8192列，原点一致；Pitch=50 μm，每Pitch 64点，物面6.4×6.4 mm。
- 保留Case1统一Anode/IJP/BM/CF膜堆和Materials；IJP全白、10 μm、n=1.5、k=0。
- D65 380–780 nm、5 nm；输出±6°/1201；沿用plane照明，保留有限源samplesPerAxis=5。
- 接收距离250 m，以满足本项目的夫琅禾费尺度；角能量效率不随此距离改变。
- input_manifest.csv记录来源/输出SHA-256、平铺倍数和材料覆盖率。
- 本轮完整仿真先使用Case1_Pitch64；本项目保留128 Pitch完整输入，运行前需评估内存。详见../../docs/README.md。
