# `CustomAperture_fraunhofer_fft.m` 使用说明

## 1. 程序简介

`CustomAperture_fraunhofer_fft.m` 是一个用于模拟自定义周期孔径夫琅禾费（Fraunhofer）远场衍射的 MATLAB 脚本。程序默认构建一个由四个圆孔组成的 `2 × 2` 晶胞，分别计算红、绿、蓝三个波长的远场衍射，并完成以下工作：

- 显示物面孔径及三个波长的远场振幅、强度分布；
- 对衍射强度进行百分位截断、对数压缩和伽马增强；
- 将三个波长映射到 RGB 三通道，生成彩色远场衍射图；
- 在启用图像照明时，以最后一个波长的强度分布作为点扩散函数（PSF），计算输入图像的衍射模糊结果；
- 统计指定观察角以内、以外的能量，并给出理论雾度。

> 本文件是脚本而不是函数。参数需要在脚本开头的“模型参量”区域直接修改。

## 2. 文件位置与依赖

主程序：

```text
src/CustomAperture_fraunhofer_fft.m
```

直接依赖的项目文件：

| 文件 | 用途 |
| --- | --- |
| `src/CustomPattern.m` | 生成旋转、平移后的圆形或椭圆形离散孔径 |
| `src/fraunhofer_fft_SPH.m` | 计算夫琅禾费衍射场 |
| `src/myFFT2_SPH.m` | 完成二维傅里叶变换及采样映射 |
| `src/interpolate.m` | 将计算结果投影/插值到目标像面 |
| `src/placeImageOnCanvas.m` | 按物理尺寸重采样输入图像并居中放置到仿真画布 |
| `src/fftconv2_same.m` | 使用 FFT 完成二维同尺寸卷积 |
| `input/USAF500.png` | `Lightflag = 1` 时的默认输入图像 |

程序会把预览图写入：

```text
ImageForShow/four_aperture_rgb_preview.png
ImageForShow/four_aperture_gray_preview.png
```

建议使用 MATLAB R2020a 或更新版本。图像读取、灰度转换、显示和缩放使用了 `rgb2gray`、`imshow`、`imresize` 等图像处理函数，通常需要 Image Processing Toolbox。

## 3. 快速运行

在 MATLAB 中将当前目录切换到项目的 `src` 文件夹，再运行脚本：

```matlab
cd('D:\03 Optical Simulation\Diffraction-Calculation\src');
CustomAperture_fraunhofer_fft
```

从 `src` 目录运行很重要，因为默认输入、输出路径均使用相对路径 `..\input` 和 `..\ImageForShow`。运行前请确认 `ImageForShow` 文件夹已经存在。

默认参数计算网格为 `2000 × 2000`，且会进行三次衍射计算、排序、绘图和卷积，可能占用较多内存并需要一定运行时间。首次试运行可先把 `nn` 从 `1000` 降为 `200` 或 `300`。

## 4. 主要参数

### 4.1 光学与采样参数

| 参数 | 默认值 | 单位 | 含义 |
| --- | ---: | --- | --- |
| `lambdalist` | `[620e-9, 550e-9, 450e-9]` | m | 三个计算波长，依次写入 RGB 的红、绿、蓝通道 |
| `z` | `0.03` | m | 传播距离 |
| `objectunit` | `1e-6` | m | 物面尺寸的便捷单位，默认对应 `1 μm` |
| `Imageunit` | `1e-3` | m | 像面尺寸的便捷单位，默认对应 `1 mm` |
| `Pitch` | `50*objectunit` | m | 单个孔径周期尺寸 |
| `r` | `10*objectunit` | m | 默认圆孔半径 |
| `nn` | `1000` | 像素/周期 | 单个周期每个方向的采样点数 |
| `Nx_period` | `2` | 周期 | 物面 x 方向周期数 |
| `Ny_period` | `2` | 周期 | 物面 y 方向周期数 |
| `Fx`, `Fy` | `20*Imageunit` | m | 像面坐标的正半轴范围；实际范围分别为 `[-Fx,Fx]`、`[-Fy,Fy]` |
| `F_size_x`, `F_size_y` | `20*Imageunit` | m | 输入图像映射到像面后的目标物理尺寸 |
| `Angle` | `2.5` | ° | 能量统计使用的目标半角 |

目标角在像面上的半径为：

```matlab
FocusR = tand(Angle) * z;
```

总网格尺寸由 `Tnn = Nx_period * nn` 决定。当前脚本使用同一个 `Tnn` 构造 x、y 两个方向的方形网格，因此建议保持 `Nx_period == Ny_period`。

### 4.2 仿真模式

#### `Lightflag`

| 值 | 模式 | 行为 |
| ---: | --- | --- |
| `0` | 平面波 | 只显示均匀输入，不执行末尾的图像卷积和雾度计算 |
| `1` | 图像照明 | 读取 `../input/USAF500.png`，执行图像卷积及雾度计算 |

如需更换输入图像，请修改：

```matlab
Input = imread('..\input\USAF500.png');
```

当前代码会调用 `rgb2gray`，因此输入文件应为 RGB 图像。若输入本身已经是灰度图，需要去掉或调整该调用。

#### `Gratingflag`

| 值 | 孔径定义方式 | 说明 |
| ---: | --- | --- |
| `0` | 同一离散单元周期复制 | 默认创建长短半轴不同的椭圆，并在 x、y 方向重复 |
| `1` | 多个离散单元拼接 | 默认构建 `2 × 2` 四孔晶胞；四个孔的旋转角依次为 `30°、0°、0°、-30°` |
| `2` | 连续匿名函数 | 定义圆孔连续函数后交给衍射函数离散化；当前实现存在兼容性限制，见“已知限制” |

当 `Gratingflag = 1` 时，`cellunit = 2`，所以 `Nx_period` 和 `Ny_period` 必须是 `2` 的整数倍，否则 `repmat` 的重复次数不是整数。默认四个孔均为圆孔，旋转圆孔不会改变结果；只有将 `a1`～`a4` 与 `b1`～`b4` 设为不同值形成椭圆时，旋转角度才会改变孔径形状。

## 5. 自定义孔径

`CustomPattern` 的调用格式为：

```matlab
Aperture = CustomPattern(angle, Xc, Yc, a, b, offsetx, offsety);
```

参数含义：

| 参数 | 含义 |
| --- | --- |
| `angle` | 椭圆逆/顺时针方向由坐标变换约定决定的旋转角，单位为度 |
| `Xc`, `Yc` | 单周期网格坐标 |
| `a`, `b` | 椭圆两个半轴；`a == b` 时为圆孔 |
| `offsetx`, `offsety` | 孔径中心相对单周期中心的偏移量，单位为 m |

例如，将第一个单元改为长半轴 `12 μm`、短半轴 `8 μm`、旋转 `30°`、x 方向偏移 `2 μm`：

```matlab
a1 = 12*objectunit;
b1 = 8*objectunit;
SinglePattern1 = CustomPattern(30, Xc, Yc, a1, b1, ...
                               2*objectunit, 0);
```

孔径必须能够被当前周期网格正确采样。一般应满足半轴和偏移后的孔径边界不超出 `[-Pitch/2, Pitch/2]`，并保证孔径直径覆盖足够多的采样点。

## 6. 计算流程

程序的主要计算顺序如下：

1. 根据 `Pitch`、周期数及 `nn` 建立物面采样网格。
2. 根据 `Fx`、`Fy` 和 `Tnn` 建立目标像面网格。
3. 根据 `Gratingflag` 构建物面振幅掩模 `Mask`。
4. 对 `lambdalist` 中的三个波长逐一调用 `fraunhofer_fft_SPH`。
5. 由复振幅 `U2` 计算振幅 `A = abs(U2)` 和强度 `I = A.^2`。
6. 对强度执行可视化增强，并依次写入 `RGB(:,:,1:3)`。
7. 显示并导出 RGB 合成图和绿色通道图。
8. 当 `Lightflag = 1` 时，用 PSF 与输入图像做 FFT 卷积，并计算目标角范围外的能量比例。

夫琅禾费衍射计算采用的基本形式为：

$$
U(X,Y) = \frac{e^{ikz}}{i\lambda z}
e^{\frac{ik}{2z}(X^2+Y^2)}
\mathcal{F}\{U_0(x,y)\}
\left(\frac{X}{\lambda z},\frac{Y}{\lambda z}\right),
$$

其中 $k=2\pi/\lambda$，`fraunhofer_fft_SPH` 使用 FFT 计算傅里叶变换项。

## 7. 强度显示增强

为了同时显示亮中心和较弱的衍射结构，程序对每个波长的强度做如下处理：

```matlab
I_norm = I ./ max(I(:));
clipLevel = 99.85% 分位强度;
I_clip = min(I_norm, clipLevel) ./ clipLevel;
I_vis = log(1 + 60*I_clip) / log(61);
I_vis = I_vis.^0.72;
```

其中：

- `0.9985` 控制高亮像素的截断百分位；
- `alpha = 60` 控制对数压缩强度；
- `gamma = 0.72` 控制暗部增强程度。

这些操作仅用于显示和 RGB 合成。物理能量与雾度计算仍使用未经该显示增强的原始强度 `I`。

## 8. 输出结果

运行后通常会产生以下窗口或文件：

| 输出 | 内容 |
| --- | --- |
| 初始振幅/孔径图 | 输入图像或平面波，以及物面孔径 `Mask` |
| 三波长结果图 | 每个波长各一幅振幅图和增强强度图 |
| RGB 合成图 | 620 nm、550 nm、450 nm 分别映射到 R、G、B |
| 绿色通道图 | `RGB(:,:,2)`，对应 550 nm 的增强强度 |
| 衍射图像 | 输入图像及其与 PSF 卷积后的结果，仅 `Lightflag = 1` 时生成 |
| 角度遮罩图 | 半径 `FocusR` 内为 1、外部为 0 |
| 命令行结果 | 目标角半径、像面范围、总能量、目标角内能量和理论雾度 |

理论雾度按目标角以外的能量占比计算：

$$
WD = 1 - \frac{I_{\mathrm{Angle}}}{I_{\mathrm{Total}}}.
$$

## 9. 常见修改示例

### 9.1 使用平面波并加快预览

```matlab
Lightflag = 0;
nn = 300;
```

### 9.2 改为单一重复椭圆孔径

```matlab
Gratingflag = 0;
a1 = 12*objectunit;
b1 = 8*objectunit;
```

注意：`a1`、`b1` 的赋值位于 `Gratingflag == 0` 分支内部，应直接修改该分支中的现有代码。

### 9.3 改变能量统计角度

```matlab
Angle = 5;
```

同时建议把命令行文本中写死的“2.5°”改成由 `Angle` 格式化输出，否则数值计算已改变，但文字标签仍显示 `2.5°`。

## 10. 仿真 D65 光源下的多波长透射衍射

### 10.1 正确的光谱叠加方法

D65 是连续光谱照明体。原程序仅计算 620 nm、550 nm 和 450 nm，并将三个结果直接放入 RGB 通道，适合快速预览，但不能代表完整的 D65 照明。

D65 各波长通常按非相干光处理，因此应叠加各波长的**强度**，不能叠加复振幅 `U2`：

$$
I_{D65}(X,Y)=
\frac{\sum_i S_{D65}(\lambda_i)T(\lambda_i)I(X,Y;\lambda_i)\Delta\lambda_i}
{\sum_i S_{D65}(\lambda_i)T(\lambda_i)\Delta\lambda_i}.
$$

其中 $S_{D65}$ 为 D65 相对光谱功率，$T$ 为样品的光谱强度透过率，$I=|U_2|^2$。若已有的是振幅透过率 $t(\lambda)$，应使用 $T(\lambda)=|t(\lambda)|^2$。

### 10.2 准备光谱数据

建议从权威 CIE 数据表取得 D65 光谱功率分布和 CIE 1931 2° 标准色度观察者颜色匹配函数，并整理为：

```text
input/D65_CIE1931_2deg.csv
```

CSV 至少包含：

| 列名 | 内容 |
| --- | --- |
| `wavelength_nm` | 波长，单位 nm |
| `D65` | D65 相对光谱功率 |
| `xbar`、`ybar`、`zbar` | CIE 1931 2° 颜色匹配函数 |

首次计算可使用 380～780 nm、10 nm 间隔，共 41 个波长；需要更高精度时再改为 5 nm 或 1 nm。

### 10.3 替换波长定义

将原来的 `lambdalist` 定义替换为：

```matlab
spectralTable = readtable('..\input\D65_CIE1931_2deg.csv');

lambdaNm = (380:10:780).';
lambdalist = lambdaNm * 1e-9;

S_D65 = interp1(spectralTable.wavelength_nm, spectralTable.D65, ...
                lambdaNm, 'pchip', 0);
xbar = interp1(spectralTable.wavelength_nm, spectralTable.xbar, ...
               lambdaNm, 'pchip', 0);
ybar = interp1(spectralTable.wavelength_nm, spectralTable.ybar, ...
               lambdaNm, 'pchip', 0);
zbar = interp1(spectralTable.wavelength_nm, spectralTable.zbar, ...
               lambdaNm, 'pchip', 0);

dLambda = gradient(lambdaNm);

% 样品强度透过率；无透射光谱时暂设为 1
Tlambda = ones(size(lambdaNm));

% 如有透射光谱，可使用：
% transmissionTable = readtable('..\input\transmission.csv');
% Tlambda = interp1(transmissionTable.wavelength_nm, ...
%                   transmissionTable.transmittance, ...
%                   lambdaNm, 'pchip', 0);
% 若透过率以百分数表示，还需执行 Tlambda = Tlambda/100;

spectralWeight = S_D65 .* Tlambda .* dLambda;
```

`transmittance` 应是 0～1 的强度透过率。光谱测量范围外采用插值值 0，可避免不合理的外推。

### 10.4 替换三波长循环

删除或注释原程序从 `figure; for c = 1:3` 开始的三波长循环，以及紧随其后的原 RGB 合成代码，改为：

```matlab
I_D65 = zeros(Tnn, Tnn);
Xcie = zeros(Tnn, Tnn);
Ycie = zeros(Tnn, Tnn);
Zcie = zeros(Tnn, Tnn);

for c = 1:numel(lambdalist)
    lambda = lambdalist(c);
    k = 2*pi/lambda;

    % 若孔径的振幅或相位响应随波长变化，应在这里更新 Mask
    [~, U2] = fraunhofer_fft_SPH( ...
        Mask, xmin, xmax, ymin, ymax, Tnn, Tnn, ...
        lambda, z, Xmin, Xmax, Ymin, Ymax, Tnn, Tnn, ...
        k, Gratingflag);

    I_lambda = abs(U2).^2;

    % 比较各波长 PSF 形状时，将单色 PSF 归一到单位能量
    I_lambda = I_lambda / max(sum(I_lambda, 'all'), eps);

    w = spectralWeight(c);
    I_D65 = I_D65 + w * I_lambda;

    % 同时累积 CIE XYZ，用于生成彩色预览
    Xcie = Xcie + w * xbar(c) * I_lambda;
    Ycie = Ycie + w * ybar(c) * I_lambda;
    Zcie = Zcie + w * zbar(c) * I_lambda;
end

I_D65 = I_D65 / max(sum(I_D65, 'all'), eps);
```

逐波长单位能量归一化适合研究宽光谱 PSF 的形状。如果需要绝对衍射效率或绝对透射功率，应去掉这一步，并统一校准入射光场、像面采样面积及探测器光谱响应。

### 10.5 显示 D65 总衍射强度

```matlab
Ishow = I_D65 / max(I_D65(:));
Ishow = log(1 + 60*Ishow) / log(61);
Ishow = Ishow.^0.72;

figure;
imagesc(X(1,:)*1e3, Y(:,1)*1e3, Ishow);
axis image;
set(gca, 'YDir', 'normal');
xlabel('x / mm');
ylabel('y / mm');
title('D65 多波长透射衍射强度');
colormap gray;
colorbar;
```

显示增强应在全部波长完成物理加权后统一执行。不能先把各波长分别归一化成显示图再相加，否则会破坏 D65 的相对光谱权重。

### 10.6 转换为 sRGB 彩色图

多波长结果不能直接塞入三个 RGB 通道。应先通过颜色匹配函数得到 XYZ，再转换为 sRGB：

```matlab
xyzScale = max(Ycie(:));
Xn = Xcie / max(xyzScale, eps);
Yn = Ycie / max(xyzScale, eps);
Zn = Zcie / max(xyzScale, eps);

% XYZ 到线性 sRGB（D65 白点）
Rlin =  3.2406*Xn - 1.5372*Yn - 0.4986*Zn;
Glin = -0.9689*Xn + 1.8758*Yn + 0.0415*Zn;
Blin =  0.0557*Xn - 0.2040*Yn + 1.0570*Zn;
RGBlin = max(cat(3, Rlin, Glin, Blin), 0);
RGBlin = RGBlin / max(max(RGBlin(:)), eps);

% sRGB 编码曲线
RGB_D65 = 12.92 * RGBlin;
idx = RGBlin > 0.0031308;
RGB_D65(idx) = 1.055 * RGBlin(idx).^(1/2.4) - 0.055;
RGB_D65 = min(max(RGB_D65, 0), 1);

figure;
imagesc(X(1,:)*1e3, Y(:,1)*1e3, RGB_D65);
axis image;
set(gca, 'YDir', 'normal');
xlabel('x / mm');
ylabel('y / mm');
title('D65 多波长透射衍射彩色结果');

exportgraphics(gcf, ...
    fullfile('..', 'ImageForShow', 'D65_multispectral_diffraction.png'), ...
    'Resolution', 200);
```

这是适合显示器观看的相对色彩预览，不等同于绝对色度测量。若要模拟具体相机，应使用相机的 R/G/B 光谱响应替代 CIE 颜色匹配函数。

### 10.7 用宽光谱 PSF 计算图像与雾度

原程序循环结束后的 `I` 只对应最后一个波长。D65 仿真中应将末尾 PSF、卷积和能量统计改为使用 `I_D65`：

```matlab
if Lightflag == 1
    PSF = I_D65 / max(sum(I_D65, 'all'), eps);
    Icanvas = placeImageOnCanvas(U0, Tnn, Tnn, ...
                                 F_size_x, F_size_y, 2*Fx, 2*Fy);
    Iout = real(fftconv2_same(Icanvas, PSF));

    figure;
    subplot(1,2,1);
    imshow(Icanvas, []);
    title('输入图像');
    subplot(1,2,2);
    imshow(Iout, []);
    title('D65 多波长衍射后图像');

    AngleMask = double(X.^2 + Y.^2 <= FocusR^2);
    I_Total = sum(I_D65, 'all');
    I_Angle = sum(I_D65 .* AngleMask, 'all');
    WD_D65 = 1 - I_Angle/I_Total;
    fprintf('D65 加权理论雾度：%g %%\n', WD_D65*100);
end
```

计算明视觉亮度加权雾度时，可在光谱权重中再乘以 `ybar`；模拟相机测量时，则应乘以对应相机通道的光谱响应。

### 10.8 精度与性能建议

- 首次调试建议使用 `nn = 200`～`400` 和 10 nm 光谱间隔；默认 `nn = 1000` 下进行 41 次衍射计算会非常耗时。
- 不要保存所有波长的完整三维强度数组；在循环中直接累加可显著降低内存占用。
- D65、透过率及颜色匹配函数必须使用一致的波长单位和覆盖范围。
- 若样品存在色散相位，必须针对每个波长更新复数 `Mask`；仅给最终强度乘 `Tlambda` 无法描述波前色散。
- 若只关心相对 PSF，可将单色 PSF 归一到单位能量；若研究绝对效率，则不要逐波长归一化，并需建立一致的辐射定标。

## 11. 已知限制与注意事项

1. **当前目录要求**：默认相对路径按 `src` 为当前目录设计；从其他目录直接执行可能找不到输入图像或把输出写到错误位置。
2. **网格默认按方阵处理**：`Tnn` 只由 `Nx_period` 计算。为避免孔径矩阵、像面网格和卷积尺寸不一致，应保持 `Nx_period == Ny_period`。
3. **四孔模式的周期约束**：`Gratingflag = 1` 时，x、y 周期数都应是 `cellunit = 2` 的整数倍。
4. **连续函数模式尚不可靠**：`fraunhofer_fft_SPH.m` 在 `flag == 2` 时调用 `discretize(Uc, xmin, ...)`，但项目实际提供的自定义函数名为 `Selfdiscretize`。若路径中没有另一个兼容的 `discretize` 实现，`Gratingflag = 2` 会报参数或类型错误。使用前应将该调用核对并改为项目提供的离散化函数。
5. **图像卷积只使用最后一个波长**：循环结束后变量 `I` 保存的是 `lambdalist(end)`，默认即 450 nm 的强度。因此 `PSF`、衍射后图像、角内能量及雾度都基于蓝光通道，不是 RGB 综合结果。
6. **显示增强不代表绝对能量**：`RGB` 各通道分别归一化、截断和增强，适合观察颜色分离，但不能直接用于比较不同波长的绝对衍射效率。
7. **像面尺寸打印符号**：代码用负值 `Xmin*2`、`Ymin*2` 输出“像面尺寸”，可能显示为负数；实际像面范围是 `[-Fx,Fx] × [-Fy,Fy]`，总宽高为 `2*Fx × 2*Fy`。
8. **局部函数 `work` 当前未启用**：其调用已被注释，且函数内部引用了未作为参数传入的 `lambda`。若取消注释调用，需要先把 `lambda` 加入函数参数或在函数内定义。
9. **内存与速度**：默认 `Tnn = 2000`，程序同时保存多个 `2000 × 2000` 双精度矩阵和 RGB 数组。增大 `nn` 或周期数会快速增加内存和计算时间。
10. **本次验证范围**：已人工核对主脚本及其直接辅助函数；本机 MATLAB R2024a 的批处理启动因 settings 插件加载失败，未完成整脚本实跑验证。

## 12. 故障排查

### 找不到 `USAF500.png`

确认 MATLAB 当前目录为项目的 `src` 文件夹，且文件存在于 `input/USAF500.png`。

### 找不到辅助函数

确保当前目录为 `src`，或将其加入 MATLAB 路径：

```matlab
addpath('D:\03 Optical Simulation\Diffraction-Calculation\src');
```

### `repmat` 提示重复次数必须为整数

在四孔模式下，把 `Nx_period` 和 `Ny_period` 设置为 `2` 的整数倍。

### 内存不足或运行过慢

先减小 `nn`，例如设为 `200`～`500`；确认结果趋势后再提高采样分辨率。

### `Gratingflag = 2` 报错

参考“已知限制”第 4 条，检查 `fraunhofer_fft_SPH.m` 中连续函数离散化所调用的函数名和参数。
