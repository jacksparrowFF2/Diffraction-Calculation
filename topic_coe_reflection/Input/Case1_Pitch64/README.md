# Case1_Pitch64：完整64×64 Pitch图片项目

- 来源：Case1六张1024×1024 PNG原样按4×4平铺，未缩放或修改像素值。
- 六张图片均为4096行×4096列；Pitch=50 μm，每Pitch 64点，物面3.2×3.2 mm。
- 统一膜堆、材料、图片原点及采样间距与Case1一致。
- D65 380–780 nm、5 nm；±6°/1201；plane照明，保留有限源samplesPerAxis=5。
- 接收距离100 m，额外高度0；input_manifest.csv记录来源及输出SHA-256。
- 本轮受资源条件选择64 Pitch做完整仿真；不将该尺寸声明为已通过128/192 Pitch的收敛要求。
- 主脚本选择projectName="Case1_Pitch64"；完整说明见../../docs/README.md。
