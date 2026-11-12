# 基础孔径衍射

本模块是项目中其它程序开展衍射计算的基础，提供孔径建模、夫琅禾费传播及相关方法参考。当前各功能模块自带私有依赖，可独立运行。

```matlab
addpath('topic_aperture_fraunhofer/src');
CustomAperture_fraunhofer_fft; % RGB 孔径衍射和原版图像卷积
% 或
CustomAperture_four_aperture_preview;
```

脚本参数在各文件开头修改。主脚本默认 `nn=1000`，首次运行可改成 200；默认周期数为 2×2。原版卷积保留既有行为，定量透射成像请使用 `topic_udc_imaging`。

- `src/private/`：孔径、夫琅禾费变换、FFT、连续函数采样、插值、图像画布与卷积，完整自带。
- `Input/USAF500.png`：默认测试图。
- `Output/`：自动创建，存放图像与预览数据，不依赖 MATLAB 当前文件夹。
- `render_four_aperture_mapping_variants.py`：读取 `Output/four_aperture_preview_data.mat`，需 Python 的 numpy/scipy/matplotlib。
- `src/archive/`：原 MATLAB 自动保存备份。

可独立复制本目录运行。图像缩放和部分显示函数需要 Image Processing Toolbox。连续函数模式现在调用自带 `Selfdiscretize`，避免误调用 MATLAB 内置同名 `discretize`。

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
