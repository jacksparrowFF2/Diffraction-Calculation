# Case2：Diamond BMP 示例图片集

这是演示几何，不代表实际产品。三套项目几何相同，仅编码格式不同。

- 主脚本设置 `projectName="Case2"`、`stack.geometry.source="images"`。
- Pitch：50 μm；periods：[16,16]；samplesPerPeriod：64。
- 图片尺寸：1024 行 × 1024 列；物面：800 × 800 μm。
- anode 白色=阳极；BM 白色=实心 BM（stack.layers.BM.invert=false）；CF_R/G/B 白色=对应膜。
- 阈值 0.5；图层文件名可不写扩展名，自动匹配本项目格式。
- 各项目 config.txt 自动提供采样、材料、观察角域和显示参数；不包含像素形状参数。
- 参考图实线 PDL 作为阳极，虚线 Mask 不作为 BM；R/B 交错、G 在四邻中心且左右交替倾斜，R:G:B 数量比 1:2:1。
- JPG 为有损格式，推荐实际二值版图优先使用 PNG/BMP；此示例使用 Quality=100。
- 材料仍取主脚本示例常数；若提供 nk CSV，请置于本项目 Materials/ 下并填写相对路径。

- IJP与其它图层使用统一图片输入：全白IJP图片，n=1.5、k=0、厚10 μm，层序Anode/IJP/BM/CF。默认801×801输出、±6°，像素内仍自动细分积分。
