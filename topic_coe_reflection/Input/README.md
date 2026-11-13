# 按项目管理输入

每个项目使用独立子文件夹，直接放置六张同尺寸、同原点的图层图片：

```text
Input/
  Case1/        # Diamond PNG 示例
    anode.png
    IJP.png      # 全白，n=1.5、k=0、厚10 μm
    BM.png
    CF_R.png
    CF_G.png
    CF_B.png
    config.txt  # 本项目的采样、材料、角域、显示配置
    README.md
  Case2/        # 同几何 BMP 示例，同样的图层名
  Case3/        # 同几何 JPG 示例，同样的图层名
  Colorimetry/  # 共享 CIE 数据，无需复制到每个项目
```

在 `src/COEReflection_main.m` 修改 `projectName="Case1"` 为其它项目名（支持中文），即可切换图片集；保持 `stack.geometry.source="images"`。

- 文件名默认 `anode`、`IJP`、`BM`、`CF_R`、`CF_G`、`CF_B`，不写扩展名会匹配 PNG/BMP/JPG/JPEG/TIF/TIFF。同一图层同时存在多种格式时需明确填写扩展名，程序不会静默选择。
- 白色表示所声明区域；示例BM白色=实心材料，invert=false；白色=开口时设stack.layers.BM.invert=true。
- 三套均为 1024×1024，Pitch=50 μm，16×16 单元，每周期 64 点；完整物面 800×800 μm。样例只用于演示，不是实测产品。
- 实际二值图样建议 PNG 或 BMP。JPG 是有损压缩，可能影响边缘，输入后会按阈值 0.5 二值化。
- 材料 CSV 可放在项目的 `Materials/` 下，主脚本的 nkFile 填写 `Materials/BM_nk.csv` 等相对路径。示例项目暂使用主脚本的占位材料参数。
- 每个项目的结果保存在 `Output/项目名称/`。改变项目并不会自动更改 Pitch、采样、材料和厚度，须按实际项目设置。
- `projectName=""` 可沿用旧的 Input/Output 根目录。

详细格式见 `../docs/README.md`。

默认自动读取各项目 config.txt，配置值覆盖主脚本的默认值。用户只导入图样，不在本程序调整像素形状。详细语法见 ../docs/README.md。

统一图层配置使用stack.layerNames和layerLevels；每项method为nk或reflectance。已提供CaseSingleNK（悬空有限膜）和CaseSingleReflectance（直接反射率）两个单图层项目，均只读取singlelayer.bmp。详见统一README的项目配置章节。

Case1~3 的默认空间光源为plane；cfg.source.type可切换point或laser。光谱可选D65、custom或monochromatic；定义与示例均见统一README。IJP图片应与其余输入同尺寸、同原点，全白表示整面有膜，透明性由k=0定义。
