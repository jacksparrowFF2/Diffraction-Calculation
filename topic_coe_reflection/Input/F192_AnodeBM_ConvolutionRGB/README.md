# F192_AnodeBM_ConvolutionRGB：三波长图样对齐

- 复制F192原始anode/BM图片，不修改原输入。原生250像素/Pitch、4Pitch单元，目标32Pitch。
- 只有Anode与BM；Anode反射率1，BM为理想不透明黑区（有效反射率0），不含IJP或CF。用于几何图样比较，不用于真实材料反射率预测。
- 620/550/450nm的RGB代理，按D65采样相对权重，非完整D65色度积分。
- 5mm均匀非相干圆盘/5cm，近轴强度卷积，wavefront="plane"。plane用于光栅远场PSF与源形状卷积；paraxial保留中心入射球面二次相位。
- 原始点源结果与卷积核保存在result.sourceConvolution，正式F192配置保持原样。详细假设见统一docs/README.md。
