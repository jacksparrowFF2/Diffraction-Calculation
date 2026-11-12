# 像素展开与批量仿真

本目录为同事编写的像素显示仿真代码，保留其像素展开和批量仿真 GUI。

```matlab
addpath('topic_pixel_unfolding/src');
OLED_PixelUnfolding_Tian_Batch_GUI;
```

保留原作者 GUI，辅助计算函数均在同一 `.m` 文件中，不依赖其他功能目录。通过界面选择待处理图片、A/B 像素模板和输出目录，输入数据由用户指定。需要 Image Processing Toolbox；Excel 报告由原程序的导出流程生成。

## 目录与独立运行

本方向统一使用以下结构：

| 目录 | 用途 |
|---|---|
| `Input/` | 输入图片序列、项目参数文本等计算数据 |
| `Output/` | 计算产生的图片、表格、MAT 等结果 |
| `src/private/` | 辅助函数；保留 MATLAB 特殊目录名，由主程序自动解析 |
| `docs/` | 使用说明、报告等文档 |
| `dist/` | 打包发布文件 |
| `src/` | 主程序源码 |

独立复制整个方向文件夹，在 MATLAB 中仅添加其 `src/` 路径后运行上述入口。输入和输出路径根据入口文件位置定位，与当前工作目录无关。辅助函数目录是命名例外：MATLAB 只有对名为 `private`、位于调用程序同一父目录下的文件夹才会自动应用私有函数解析规则。因此辅助函数放在 `src/private/`，保留小写名称，不改为 `PrivateFunction`，也不单独 `addpath`。
