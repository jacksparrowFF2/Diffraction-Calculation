# F192 Anode/BM 平面波结构对比

anode.png、BM.png 原样复制自 F192_AnodeBM_ConvolutionRGB，无缩放、裁切或轮廓修改。两图均为1000×1000，对应4×4 Pitch，Pitch=59.06 μm、每Pitch 250点。

仅使用Anode/BM：Anode强度反射率1、相位180°，BM反相后作为上层零反射不透明区域；不含CF/IJP。平面波，620/550/450 nm，D65三点相对权重，仅为RGB代理。

主入口 COEReflection_main.m 选择本项目时默认开启 runPlaneWaveComparison：完整4×4 Pitch局部图片与16×16 Pitch有限阵列（各轴4次重复）对比，统一±6°角域、1201×1201网格、共同曝光。局部图片是完整输入单元，不声称为最小晶胞。输出在 Output/F192_AnodeBM_PlaneRGB/PlaneWaveComparison。

完整使用说明见 ../../docs/README.md 的“真实Anode/BM平面波局部单元与阵列对比”。
