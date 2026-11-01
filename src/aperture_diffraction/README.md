# 基础孔径衍射

```matlab
addpath('src/aperture_diffraction');
CustomAperture_fraunhofer_fft; % RGB 孔径衍射和原版图像卷积
% 或
CustomAperture_four_aperture_preview;
```

脚本参数在各文件开头修改。主脚本默认 `nn=1000`，首次运行可改成 200；默认周期数为 2×2。原版卷积保留既有行为，定量透射成像请使用 `udc_imaging`。

- `private/`：孔径、夫琅禾费变换、FFT、连续函数采样、插值、图像画布与卷积，完整自带。
- `input/USAF500.png`：默认测试图。
- `output/`：自动创建，存放图像与预览数据，不依赖 MATLAB 当前文件夹。
- `render_four_aperture_mapping_variants.py`：读取 `output/four_aperture_preview_data.mat`，需 Python 的 numpy/scipy/matplotlib。
- `archive/`：原 MATLAB 自动保存备份。

可独立复制本目录运行。图像缩放和部分显示函数需要 Image Processing Toolbox。连续函数模式现在调用自带 `Selfdiscretize`，避免误调用 MATLAB 内置同名 `discretize`。
