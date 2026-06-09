%% 验证叠加非周期相位去相关层后的夫琅禾费衍射
clear, clc, close all
%% 模型参量
% 定义波长，e-9代表纳米
lambdalist = [620e-9,550e-9,450e-9];
% 传播距离,模拟人眼视距2~40cm
z = 0.2;
% 设定物面单元尺寸
objectunit = 1e-6;
% 设定像面单元尺寸
Imageunit = 1e-3;
% 设定周期尺寸,直观单位微米
Pitch = 50*objectunit;
% 圆孔半径,e-6代表μm
r = 10*objectunit;
% 定义物面的光栅类型
% flag=0代表使用相同离散矩阵（自定义的形状或者图片）
% flag=1代表使用不相同的离散矩阵
% flag=2代表使用连续函数
Gratingflag = 1;
% 定义是否仿真像面的衍射效果，平面波为0，图片为1
Lightflag = 0;
% 分辨率，设定单个周期的分辨率
nn = 1000;
% 设定物面X、Y方向周期
Nx_period = 2;         
Ny_period = 2;         
% 设定目标像面单边尺寸，直观单位mm
Fx =40*Imageunit;
Fy =40*Imageunit;
% 设定物体在像面的尺寸，直观单位mm
F_size_x = 20*Imageunit;   % F在像面上的物理宽度
F_size_y = 20*Imageunit;   % F在像面上的物理高度
%设定目标角，用于寻找指定角度之外和之内的强度和，单位是°，黄度是2.5°
Angle = 2.5;
% 计算目标角对应传播距离下的长度
FocusR = tand(Angle)*z;
fprintf('像面上2.5°范围半径为：%g mm\n', FocusR*1e3);

% 设定单个周期的最大与最小值
x_cell_min = -Pitch/2;
x_cell_max = Pitch/2;
y_cell_min = -Pitch/2;
y_cell_max = Pitch/2;
% 设定单个周期的格点坐标
x_cell = linspace(x_cell_min,x_cell_max,nn);
y_cell = linspace(y_cell_min,y_cell_max,nn);

[Xc, Yc] = meshgrid(x_cell, y_cell);


%% 计算物面及像面参数
% 物面尺寸x方向总最小值
xmin = -Nx_period*Pitch/2;
% 物面尺寸x方向总最大值
xmax = -xmin;
% 物面尺寸y方向总最小值
ymin = -Ny_period*Pitch/2;
% 物面尺寸y方向总最大值
ymax = -ymin;

% 像面与物面尺寸之比
ImageMaskRatioX = Fx/xmax;
ImageMaskRatioY = Fy/ymax;

% 像面尺寸x方向最小值
Xmin = xmin*ImageMaskRatioX;
% 像面尺寸x方向最大值
Xmax = -Xmin;
% 像面尺寸y方向最小值
Ymin = ymin*ImageMaskRatioY;
% 像面尺寸y方向最大值
Ymax = -Ymin;

% 输出像面尺寸
fprintf('像面尺寸为：%g x %g mm\n', Xmin*2*1e3, Ymin*2*1e3);
% 计算总的分辨率
Tnn = Nx_period*nn;
% 
% 构建物面x方向矢量
x = linspace(xmin, xmax, Tnn);
% 构建物面y方向矢量
y = linspace(ymin, ymax, Tnn);
% 构建物面xy离散坐标点
[x, y] = meshgrid(x, y);
% 构建像面X方向矢量
X = linspace(Xmin, Xmax, Tnn);
% 构建像面Y方向矢量
Y = linspace(Ymin, Ymax, Tnn);
% 构建像面XY离散坐标点
[X, Y] = meshgrid(X, Y);


% 给结果矩阵预分配内存
RGB = zeros(Tnn, Tnn, 3);
RGB_Origin = zeros(Tnn, Tnn, 3);
RGB_NoPhase = zeros(Tnn, Tnn, 3);
GreenCompare_NoPhase = [];
GreenCompare_WithPhase = [];
GreenIntensity_NoPhase = [];
GreenIntensity_WithPhase = [];
PeakIntensity_NoPhase = zeros(1, numel(lambdalist));
PeakIntensity_WithPhase = zeros(1, numel(lambdalist));
ZeroOrderEnergy_NoPhase = zeros(1, numel(lambdalist));
ZeroOrderEnergy_WithPhase = zeros(1, numel(lambdalist));
ZeroOrderRetention = zeros(1, numel(lambdalist));

%% 非周期相位去相关层参数
% 本脚本模拟 COE OLED 息屏/黑态下的环境反射光，因此只考虑双程相位：
% phi = 4*pi/lambda*(n_high - n_low)*h(x,y)
n_high = 1.75;
n_low = 1.55;
delta_n = n_high - n_low;
phaseHeight = 220e-9;          % 高折圆形微岛高度，推荐从 200~250 nm 起步
zeroOrderRadiusFactor = 0.45;  % 0 级主斑判读半径，按各波长一阶衍射间距的比例估算
useBestDoeDesign = true;       % true 时自动读取 DOE 最优候选；读取失败则使用下方手动参数

% 以下参数作为手动兜底值；自动导入成功时会被 DOE 最优候选覆盖。
phaseIslandRadiusList = [1.0, 1.5] * 1e-6; % 多半径圆形高折微岛半径，单位 m
phaseFillFactor = 0.40;                    % 目标填充率，实际填充率由离散采样和拒绝采样共同决定
phaseMinDistanceFactor = 1.00;             % 最小中心距系数，乘以 2*最大半径
phaseMinDistance = phaseMinDistanceFactor * 2 * max(phaseIslandRadiusList);
phaseSeed = 5;                             % DOE 最优候选对应 seed

if useBestDoeDesign
    try
        bestDoeDesign = loadBestMultiRadiusPhaseDesign();
        if isfield(bestDoeDesign, 'phaseHeight')
            phaseHeight = bestDoeDesign.phaseHeight;
        end
        phaseIslandRadiusList = bestDoeDesign.phaseIslandRadiusList;
        phaseFillFactor = bestDoeDesign.phaseFillFactor;
        phaseMinDistanceFactor = bestDoeDesign.phaseMinDistanceFactor;
        phaseMinDistance = bestDoeDesign.phaseMinDistance;
        phaseSeed = bestDoeDesign.phaseSeed;
        fprintf('已导入 DOE 最优设计：%s\n', bestDoeDesign.sourcePath);
        fprintf('DOE phaseHeight = %.0f nm\n', phaseHeight * 1e9);
        fprintf('DOE Objective = %.4g\n', bestDoeDesign.objective);
    catch ME
        warning('DOEImport:LoadFailed', 'DOE 最优设计导入失败，使用脚本内手动参数。原因：%s', ME.message);
    end
end

if Lightflag ==1
    % 导入要被衍射的图片
    Input = imread('..\input\USAF500.png');
    U0 = double(rgb2gray(Input));
    figure(1);hold on;
    subplot(1,2,1);
    imshow(U0);
    title('初始振幅'); 
    axis on;axis image;xlabel('x');ylabel('y')
else
    U0 = 1; %不参与最终计算，输入光是平面波
    subplot(1,2,1);
    imshow(U0*ones(Nx_period*nn,Ny_period*nn));
    title('初始振幅'); 
    axis on;axis image;xlabel('x');ylabel('y')
end


if Gratingflag ==0
    % 创建离散的单元函数，长半轴为a，短半轴为b
    a1 = r/0.8; b1=r;
    SinglePattern = CustomPattern(0,Xc,Yc,a1,b1,0,0);
    % 在x和y方向上分别重复Nx_period、Ny_period次周期
    Mask = repmat(SinglePattern, Ny_period, Nx_period);
    % figure;
    subplot(1,2,2);
    imshow(Mask, []);
    title('光栅单个Pitch内形状: Circular Aperture'); 
    axis on;axis image;xlabel('x');ylabel('y')
elseif Gratingflag ==1
    % 晶胞
    cellunit = 2;
    % 创建离散的单元函数1，长半轴为a，短半轴为b
    a1 = r; b1=r;
    SinglePattern1 = CustomPattern(30,Xc,Yc,a1,b1,0,0);
    % 创建离散的单元函数2，长半轴为a，短半轴为b
    a2 = r; b2=r;
    SinglePattern2 = CustomPattern(0,Xc,Yc,a2,b2,0,0);
    % 创建离散的单元函数3，长半轴为a，短半轴为b
    a3 = r; b3=r;
    SinglePattern3 = CustomPattern(0,Xc,Yc,a3,b3,0,0);
    % 创建离散的单元函数4，长半轴为a，短半轴为b
    a4 = r; b4=r;
    SinglePattern4 = CustomPattern(-30,Xc,Yc,a4,b4,0,0);
    % 拼接晶胞矩阵
    cell = [SinglePattern1 SinglePattern2;SinglePattern3 SinglePattern4];
    Mask = repmat(cell, Ny_period/cellunit, Nx_period/cellunit);
    subplot(1,2,2);
    imshow(Mask, []);
    title('光栅单个Pitch内形状: Circular Aperture'); 
    axis on;axis image;xlabel('x');ylabel('y')
elseif Gratingflag ==2
    % 创建用于作为物面的函数，一定是连续函数
    Mask = @(Xc, Yc) Xc.^2 + Yc.^2 <= r^2;
end

% 相位层需要和物面孔径在同一张离散网格上相乘。
% 如果原孔径是连续函数，这里先离散化；后续 FFT 统一按离散复振幅处理。
if Gratingflag == 2
    BaseMask = discretize(Mask, xmin, xmax, ymin, ymax, Tnn, Tnn);
else
    BaseMask = Mask;
end
FFTflag = 0;

% 生成非周期泊松盘分布的多半径圆形高折微岛高度图。
% h_map 的单位是 m，只有微岛区域为 phaseHeight，其余区域为 0。
[h_map, phaseCenters, actualFillFactor] = generateMultiRadiusPoissonDiskHeightMap( ...
    x, y, xmin, xmax, ymin, ymax, phaseHeight, phaseIslandRadiusList, ...
    phaseMinDistance, phaseFillFactor, phaseSeed);
fprintf('相位层圆形微岛数量：%d\n', size(phaseCenters, 1));
fprintf('相位层圆形微岛半径集合：%s um\n', mat2str(phaseIslandRadiusList * 1e6));
fprintf('相位层实际填充率：%.2f %%\n', actualFillFactor * 100);

% 用 550 nm 作为代表波长显示双程相位 map；实际计算时 RGB 三个波长会分别换算相位。
lambda_phase_show = 550e-9;
phi_double_show = 4 * pi / lambda_phase_show * delta_n * h_map;

figure;
subplot(1, 2, 1);
imagesc(x(1,:) * 1e6, y(:,1) * 1e6, h_map * 1e9);
axis image;
set(gca, 'YDir', 'normal');
xlabel('x / um'); ylabel('y / um');
title('非周期泊松盘圆形高折微岛高度图 / nm');
colormap gray; colorbar;

subplot(1, 2, 2);
scatter(phaseCenters(:,1) * 1e6, phaseCenters(:,2) * 1e6, ...
    12 + 20 * phaseCenters(:,3) / max(phaseCenters(:,3)), phaseCenters(:,3) * 1e6, 'filled');
axis image;
xlim([xmin xmax] * 1e6); ylim([ymin ymax] * 1e6);
xlabel('x / um'); ylabel('y / um');
title('多半径圆形微岛中心泊松盘分布');
grid on;
colorbar;

figure;
subplot(1, 2, 1);
imagesc(x(1,:) * 1e6, y(:,1) * 1e6, phi_double_show / pi);
axis image;
set(gca, 'YDir', 'normal');
xlabel('x / um'); ylabel('y / um');
title('550 nm 双程相位 map / pi');
colormap gray; colorbar;

subplot(1, 2, 2);
surf(x * 1e6, y * 1e6, phi_double_show / pi, 'EdgeColor', 'none');
axis tight;
xlabel('x / um'); ylabel('y / um'); zlabel('\phi / \pi');
title('550 nm 双程相位 3D 分布');
colormap gray; colorbar;
view(42, 32);
exportgraphics(gcf, fullfile('..', 'ImageForShow', 'aperiodic_phase_map_2d_3d.png'), 'Resolution', 200);

% 连续函数形式
% circle2 = @(x, y) (x-10e-6).^2 + (y-10e-6).^2 <= r^2;
% circle = @(x, y)circle1(x,y)&circle2(x,y);
figure;
for c = 1:3
    lambda = lambdalist(c);
    % 计算波矢量
    k = 2 * pi / lambda;
    % 先计算无相位层的基准远场，便于和叠加相位层后的结果做同尺度对比。
    [~, U2_NoPhase] = fraunhofer_fft_SPH(BaseMask, xmin, xmax, ymin, ymax, Tnn, Tnn, ...
                            lambda, z, Xmin, Xmax, Ymin, Ymax, Tnn, Tnn, ...
                            k,FFTflag);
    % 反射光双程穿过相位层，故相位延迟是单程透射的 2 倍。
    phi_double = 4 * pi / lambda * delta_n * h_map;
    PhaseMask = exp(1j * phi_double);
    ComplexMask = BaseMask .* PhaseMask;
    % 数值解-优先计算数值解
    [~, U2_WithPhase] = fraunhofer_fft_SPH(ComplexMask, xmin, xmax, ymin, ymax, Tnn, Tnn, ...
                            lambda, z, Xmin, Xmax, Ymin, Ymax, Tnn, Tnn, ...
                            k,FFTflag);
    % 绘制图像
    % 计算振幅
    A = abs(U2_WithPhase);
    % 强度分布
    I =A.^2;
    I_NoPhase = abs(U2_NoPhase).^2;
    PeakIntensity_NoPhase(c) = max(I_NoPhase(:));
    PeakIntensity_WithPhase(c) = max(I(:));
    zeroOrderRadius = zeroOrderRadiusFactor * lambda * z / Pitch;
    zeroOrderMask = X.^2 + Y.^2 <= zeroOrderRadius^2;
    ZeroOrderEnergy_NoPhase(c) = sum(I_NoPhase(zeroOrderMask), 'all');
    ZeroOrderEnergy_WithPhase(c) = sum(I(zeroOrderMask), 'all');
    ZeroOrderRetention(c) = ZeroOrderEnergy_WithPhase(c) / ZeroOrderEnergy_NoPhase(c);
    % % 对数压缩亮部
    % alpha = 100;
    % I_log = log(1 + alpha * I);
    % I_log = I_log ./ max(I_log(:));
    % % 伽马增强暗部
    % gamma = 1;
    % I_enhanced = I_log .^ gamma;

    
    I_norm = I ./ max(I(:));
    I_sorted = sort(I_norm(:));
    clipLevel = I_sorted(max(1, round(0.9985 * numel(I_sorted))));
    I_clip = min(I_norm, clipLevel) ./ clipLevel;
    alpha = 60;
    gamma = 0.72;
    I_vis = log(1 + alpha * I_clip) / log(1 + alpha);
    I_vis = I_vis .^ gamma;

    subplot(3, 2, 2*c-1)
    imagesc(X(1,:), Y(:,1), A)
    xlabel("x");ylabel("y");zlabel("|U|")
    axis image;title("光栅振幅图像"+num2str(lambda*1e9) +'nm')
    colormap gray;
    
    subplot(3, 2, 2*c)
    imagesc(X(1,:), Y(:,1), I_vis)
    xlabel("x");ylabel("y");zlabel("|U|")
    axis image;title("光栅强度分布"+num2str(lambda*1e9) +'nm')
    colormap gray;
    
    % RGB(:,:,c) = (I./max(max(I))).^gamma;
    % 原版对比度
    RGB_Origin(:,:,c) = A;
    % 对比度增强版
    RGB(:,:,c) = I_vis;

    I_norm_no_phase = I_NoPhase ./ max(I_NoPhase(:));
    I_sorted_no_phase = sort(I_norm_no_phase(:));
    clipLevel_no_phase = I_sorted_no_phase(max(1, round(0.9985 * numel(I_sorted_no_phase))));
    I_clip_no_phase = min(I_norm_no_phase, clipLevel_no_phase) ./ clipLevel_no_phase;
    I_vis_no_phase = log(1 + alpha * I_clip_no_phase) / log(1 + alpha);
    I_vis_no_phase = I_vis_no_phase .^ gamma;
    RGB_NoPhase(:,:,c) = I_vis_no_phase;

    if c == 2
        % 保存绿色通道用于相位增加前后的直观对比。
        commonMax = max([I_NoPhase(:); I(:)]);
        GreenCompare_NoPhase = enhanceIntensityForDisplay(I_NoPhase / commonMax, alpha, gamma);
        GreenCompare_WithPhase = enhanceIntensityForDisplay(I / commonMax, alpha, gamma);
        GreenIntensity_NoPhase = I_NoPhase;
        GreenIntensity_WithPhase = I;
    end
    
end

fprintf('\n相位增加前后峰值强度对比：\n');
for c = 1:numel(lambdalist)
    fprintf('%g nm: 无相位 %.4g, 有相位 %.4g, 峰值比例 %.2f %%，0级主斑保持 %.2f %%\n', ...
        lambdalist(c) * 1e9, PeakIntensity_NoPhase(c), PeakIntensity_WithPhase(c), ...
        PeakIntensity_WithPhase(c) / PeakIntensity_NoPhase(c) * 100, ...
        ZeroOrderRetention(c) * 100);
end

%% RGB合成显示色分离
figure;
subplot(1, 2, 1);
imagesc(X(1,:)*1e3, Y(:,1)*1e3, RGB_NoPhase);
axis image;
xlabel('x / mm');
ylabel('y / mm');
title('RGB Far-field Without Phase Layer');
set(gca, 'YDir', 'normal');
clim([0 1]);

subplot(1, 2, 2);
imagesc(X(1,:)*1e3, Y(:,1)*1e3, RGB);
axis image;
xlabel('x / mm');
ylabel('y / mm');
title('RGB Far-field With Aperiodic Phase Layer');
set(gca, 'YDir', 'normal');
clim([0 1]);
exportgraphics(gcf, fullfile('..', 'ImageForShow', 'four_aperture_phase_before_after_rgb.png'), 'Resolution', 200);

figure;
subplot(1, 2, 1);
imagesc(X(1,:)*1e3, Y(:,1)*1e3, GreenCompare_NoPhase);
axis image;
xlabel('x / mm');
ylabel('y / mm');
title('Green Channel Without Phase Layer');
set(gca, 'YDir', 'normal');
colormap gray;
colorbar;
clim([0 1]);

subplot(1, 2, 2);
imagesc(X(1,:)*1e3, Y(:,1)*1e3, GreenCompare_WithPhase);
axis image;
xlabel('x / mm');
ylabel('y / mm');
title('Green Channel With Aperiodic Phase Layer');
set(gca, 'YDir', 'normal');
colormap gray;
colorbar;
clim([0 1]);
exportgraphics(gcf, fullfile('..', 'ImageForShow', 'four_aperture_phase_before_after_green.png'), 'Resolution', 200);

% ① Figure 7：绿色通道相位增加前后 dB 对比。
% 使用同一个参考最大值和同一个色标，用于观察非零级强衍射峰是否被削弱，
% 同时检查 0 级主斑是否仍集中在中心标记圈内。
greenReference = max(GreenIntensity_NoPhase(:));
greenDbFloor = -60;
greenLambda = lambdalist(2);
greenZeroOrderRadius = zeroOrderRadiusFactor * greenLambda * z / Pitch;
greenZeroOrderMask = X.^2 + Y.^2 <= greenZeroOrderRadius^2;
greenZeroOrderRetention = sum(GreenIntensity_WithPhase(greenZeroOrderMask), 'all') / ...
    sum(GreenIntensity_NoPhase(greenZeroOrderMask), 'all');
GreenDb_NoPhase = intensityToDb(GreenIntensity_NoPhase, greenReference, greenDbFloor);
GreenDb_WithPhase = intensityToDb(GreenIntensity_WithPhase, greenReference, greenDbFloor);

figure;
subplot(1, 2, 1);
imagesc(X(1,:)*1e3, Y(:,1)*1e3, GreenDb_NoPhase);
axis image;
xlabel('x / mm');
ylabel('y / mm');
title('Green Without Phase / dB');
set(gca, 'YDir', 'normal');
colormap turbo;
colorbar;
clim([greenDbFloor 0]);
hold on;
plotCenterCircle(greenZeroOrderRadius * 1e3, 'w--');

subplot(1, 2, 2);
imagesc(X(1,:)*1e3, Y(:,1)*1e3, GreenDb_WithPhase);
axis image;
xlabel('x / mm');
ylabel('y / mm');
title(sprintf('Green With Phase / dB, 0-order %.1f%%', greenZeroOrderRetention * 100));
set(gca, 'YDir', 'normal');
colormap turbo;
colorbar;
clim([greenDbFloor 0]);
hold on;
plotCenterCircle(greenZeroOrderRadius * 1e3, 'w--');
exportgraphics(gcf, fullfile('..', 'ImageForShow', 'green_before_after_unified_db.png'), 'Resolution', 200);

% ② Figure 8：绿色通道差分图。
% 正值表示相位层增强该角度能量，负值表示削弱；需区分能量是回到 0 级主斑，
% 还是被推向非零级峰周围或宽角背景。
GreenDiff = (GreenIntensity_WithPhase - GreenIntensity_NoPhase) / greenReference;
diffAbs = sort(abs(GreenDiff(:)));
diffLimit = diffAbs(max(1, round(0.995 * numel(diffAbs))));
diffLimit = max(diffLimit, eps);

figure;
imagesc(X(1,:)*1e3, Y(:,1)*1e3, GreenDiff);
axis image;
xlabel('x / mm');
ylabel('y / mm');
title(sprintf('Green Difference, 0-order retention %.1f%%', greenZeroOrderRetention * 100));
set(gca, 'YDir', 'normal');
colormap(makeRedBlueCmap(256));
colorbar;
clim([-diffLimit diffLimit]);
hold on;
plotCenterCircle(greenZeroOrderRadius * 1e3, 'k--');
exportgraphics(gcf, fullfile('..', 'ImageForShow', 'green_phase_difference.png'), 'Resolution', 200);

% ③ Figure 9：中心水平截面强度曲线。
% 上图线性归一化，下图 dB；用于定量比较中心水平线上的非零级峰值削弱、
% 0 级主斑保持和旁瓣/背景变化。
centerRow = round(size(GreenIntensity_NoPhase, 1) / 2);
xLine_mm = X(centerRow, :) * 1e3;
lineNoPhase = GreenIntensity_NoPhase(centerRow, :) / greenReference;
lineWithPhase = GreenIntensity_WithPhase(centerRow, :) / greenReference;
lineNoPhaseDb = intensityToDb(GreenIntensity_NoPhase(centerRow, :), greenReference, greenDbFloor);
lineWithPhaseDb = intensityToDb(GreenIntensity_WithPhase(centerRow, :), greenReference, greenDbFloor);

figure;
subplot(2, 1, 1);
plot(xLine_mm, lineNoPhase, 'k-', 'LineWidth', 1.2); hold on;
plot(xLine_mm, lineWithPhase, 'r-', 'LineWidth', 1.2);
grid on;
xlabel('x / mm');
ylabel('Normalized intensity');
title('Green Center Horizontal Section / Linear');
legend('Without phase', 'With phase', 'Location', 'best');
xline(-greenZeroOrderRadius * 1e3, 'b--', '0-order edge');
xline(greenZeroOrderRadius * 1e3, 'b--', '0-order edge');

subplot(2, 1, 2);
plot(xLine_mm, lineNoPhaseDb, 'k-', 'LineWidth', 1.2); hold on;
plot(xLine_mm, lineWithPhaseDb, 'r-', 'LineWidth', 1.2);
grid on;
xlabel('x / mm');
ylabel('Intensity / dB');
title(sprintf('Green Center Horizontal Section / dB, 0-order %.1f%%', greenZeroOrderRetention * 100));
ylim([greenDbFloor 0]);
legend('Without phase', 'With phase', 'Location', 'best');
xline(-greenZeroOrderRadius * 1e3, 'b--', '0-order edge');
xline(greenZeroOrderRadius * 1e3, 'b--', '0-order edge');
exportgraphics(gcf, fullfile('..', 'ImageForShow', 'green_center_horizontal_section_linear_db.png'), 'Resolution', 200);

figure;
imagesc(X(1,:)*1e3, Y(:,1)*1e3, RGB);
axis image;
xlabel('x / mm');
ylabel('y / mm');
title('RGB Far-field Color Separation With Aperiodic Phase Layer');
set(gca, 'YDir', 'normal');
clim([0 1]);
exportgraphics(gcf, fullfile('..', 'ImageForShow', 'four_aperture_with_phase_rgb_preview.png'), 'Resolution', 200);


figure;
imagesc(X(1,:)*1e3, Y(:,1)*1e3, RGB(:,:,2));
axis image;
xlabel('x / mm');
ylabel('y / mm');
title('Four-aperture diffraction with aperiodic phase layer, green channel');
set(gca, 'YDir', 'normal');
colormap gray;
colorbar;
clim([0 1]);
exportgraphics(gcf, fullfile('..', 'ImageForShow', 'four_aperture_with_phase_gray_preview.png'), 'Resolution', 200);

% work(X, Y, U2, "数值解");

%% 计算图像衍射结果
if Lightflag ==1
    % 计算点扩散函数
    PSF = I;
    % 归一化
    PSF = PSF / max(PSF(:));
    % 创建与PSF相同的画布，输入图片重采样
    Icanvas =placeImageOnCanvas(U0, Tnn, Tnn, F_size_x, F_size_y, 2*Fx, 2*Fy);
    
    figure;
    subplot(1,2,1);
    imshow(Icanvas,[]);
    title('输入图像'); 
    axis on;axis image;xlabel('x');ylabel('y')
    Iout = fftconv2_same(Icanvas, PSF);

    subplot(1,2,2);
    imshow(Iout,[]);
    title('衍射后图像'); 
    axis on;axis image;xlabel('x');ylabel('y')


    % 雾度计算
    % 构建像面上在指定角度下的范围遮罩
    AngleMask = double(X.^2 + Y.^2 <= FocusR^2);
    figure;
    imagesc(X(1,:),Y(:,1),AngleMask)
    title('像面角度遮罩: Circular Aperture'); 
    axis on;axis image;xlabel('x');ylabel('y')
    
    %计算总能量（因为是计算比例，所以没有乘以像面xy步长）
    I_Total = sum(I,'all');
    fprintf('总能量：%g\n', I_Total);
    % 计算2.5度范围内的能量
    I_Angle = sum(I.*AngleMask,'all');
    fprintf('2.5度范围内总能量：%g\n', I_Angle);
    %计算雾度
    WD = 1-I_Angle/I_Total;
    fprintf('理论雾度：%g %%\n', WD*100);
end 




%% 局部函数
function [hMap, centers, actualFillFactor] = generateMultiRadiusPoissonDiskHeightMap( ...
    xGrid, yGrid, xmin, xmax, ymin, ymax, islandHeight, islandRadius, ...
    minCenterDistance, targetFillFactor, randomSeed)
    %GENERATEMULTIRADIUSPOISSONDISKHEIGHTMAP 生成非周期泊松盘多半径圆形微岛高度图
    %   这里的微岛是可量产友好的圆形 primitive：
    %   1. 每个微岛俯视图为圆形；
    %   2. 中心点满足最小距离约束，避免形成规则周期；
    %   3. 输出 hMap 后可直接换算为双程相位 phi(x,y)。

    rng(randomSeed);
    radiusList = islandRadius(:)';
    maxRadius = max(radiusList);
    xVector = xGrid(1, :);
    yVector = yGrid(:, 1);
    hMap = zeros(size(xGrid));

    validXmin = xmin + maxRadius;
    validXmax = xmax - maxRadius;
    validYmin = ymin + maxRadius;
    validYmax = ymax - maxRadius;
    if validXmin >= validXmax || validYmin >= validYmax
        error('圆形微岛最大半径过大，已经超过当前物面尺寸。');
    end

    domainArea = (xmax - xmin) * (ymax - ymin);
    meanIslandArea = mean(pi * radiusList.^2);
    targetIslandCount = ceil(targetFillFactor * domainArea / meanIslandArea);
    maxAttempts = max(20000, targetIslandCount * 1000);

    centers = zeros(targetIslandCount, 3);
    centerCount = 0;
    coveredArea = 0;
    targetArea = targetFillFactor * domainArea;
    attemptCount = 0;

    % 简单拒绝采样版泊松盘。当前 100 um 级窗口、几百个微岛时足够直观稳定。
    while coveredArea < targetArea && attemptCount < maxAttempts
        attemptCount = attemptCount + 1;
        candidateRadius = radiusList(randi(numel(radiusList)));
        candidateX = validXmin + rand() * (validXmax - validXmin);
        candidateY = validYmin + rand() * (validYmax - validYmin);

        if centerCount == 0
            acceptCandidate = true;
        else
            dx = centers(1:centerCount, 1) - candidateX;
            dy = centers(1:centerCount, 2) - candidateY;
            acceptCandidate = all(dx.^2 + dy.^2 >= minCenterDistance^2);
        end

        if acceptCandidate
            centerCount = centerCount + 1;
            if centerCount > size(centers, 1)
                centers(end + targetIslandCount, :) = 0; %#ok<AGROW>
            end
            centers(centerCount, :) = [candidateX, candidateY, candidateRadius];
            coveredArea = coveredArea + pi * candidateRadius^2;
        end
    end

    centers = centers(1:centerCount, :);
    if coveredArea < targetArea
        fprintf('警告：泊松盘采样未达到目标填充率，目标面积 %.4g，实际面积 %.4g。\n', ...
            targetArea, coveredArea);
    end

    % 将矢量化圆形微岛栅格化成高度图；只更新局部包围盒，避免全图逐岛扫描。
    for idx = 1:centerCount
        cx = centers(idx, 1);
        cy = centers(idx, 2);
        islandRadius = centers(idx, 3);

        xIndex = find(abs(xVector - cx) <= islandRadius);
        yIndex = find(abs(yVector - cy) <= islandRadius);
        [localX, localY] = meshgrid(xVector(xIndex), yVector(yIndex));
        localCircle = (localX - cx).^2 + (localY - cy).^2 <= islandRadius^2;

        localHeight = hMap(yIndex, xIndex);
        localHeight(localCircle) = islandHeight;
        hMap(yIndex, xIndex) = localHeight;
    end

    actualFillFactor = nnz(hMap > 0) / numel(hMap);
end

function I_vis = enhanceIntensityForDisplay(I_norm, alpha, gamma)
    %ENHANCEINTENSITYFORDISPLAY 对强度图做统一的对数压缩和伽马增强
    %   输入 I_norm 应已经按共同参考强度归一化，便于相位增加前后做同尺度对比。
    I_norm = max(I_norm, 0);
    I_clip = min(I_norm, 1);
    I_vis = log(1 + alpha * I_clip) / log(1 + alpha);
    I_vis = I_vis .^ gamma;
end

function I_dB = intensityToDb(I, referenceIntensity, dBFloor)
    %INTENSITYTODB 将强度转换为 dB，并设置下限，便于统一色标显示。
    I_safe = max(I, referenceIntensity * 10^(dBFloor / 10));
    I_dB = 10 * log10(I_safe / referenceIntensity);
end

function cmap = makeRedBlueCmap(n)
    %MAKEREDBLUECMAP 生成红-白-蓝发散色图，用于显示正负差分。
    if nargin < 1
        n = 256;
    end

    halfN = floor(n / 2);
    bluePart = [linspace(0, 1, halfN)', linspace(0, 1, halfN)', ones(halfN, 1)];
    redPart = [ones(n - halfN, 1), linspace(1, 0, n - halfN)', linspace(1, 0, n - halfN)'];
    cmap = [bluePart; redPart];
end

function plotCenterCircle(radius_mm, lineSpec)
    %PLOTCENTERCIRCLE 在当前坐标轴上标记 0 级主斑判读半径。
    theta = linspace(0, 2*pi, 360);
    plot(radius_mm * cos(theta), radius_mm * sin(theta), lineSpec, 'LineWidth', 1.0);
end

function work(x, y, U, name)
    % 振幅
    A = abs(U);
    % 计算强度和点扩散函数
    I =A.^2;
    % 归一化
    I = I/max(I(:));
    
    % 对数压缩亮部
    alpha = 100;
    I_log = log(1 + alpha * I);
    I_log = I_log ./ max(I_log(:));
    
    % 伽马增强暗部
    gamma = 0.5;
    I_enhanced = I_log .^ gamma;
    
    figure
    subplot(1, 2, 1)
    % surf(x, y, abs(U), 'EdgeColor', 'none', 'FaceAlpha', 0.8)
    imagesc(x(1,:), y(:,1), A)
    xlabel("x")
    ylabel("y")
    zlabel("|U|")
    % axis equal
    axis image
    title(name+"光栅振幅图像"+num2str(lambda*1e9) +'nm')
    colormap gray;
    cmin = min(A(:)); % 数据最小值
    cmax = max(A(:)); % 数据最大值
    colorbar('Limits', [cmin, cmax]); % 固定颜色范围
    
    
    subplot(1, 2, 2)
    % subplot(2, 2, 4)
    imagesc(x(1,:), y(:,1), I_enhanced)
    xlabel("x")
    ylabel("y")
    zlabel("|U|^2")
    axis image
    % axis equal
    title(name+"光栅强度分布"+num2str(lambda*1e9) +'nm')
    colormap gray;
    cmin = min(I(:)); % 数据最小值
    cmax = max(I(:)); % 数据最大值
    colorbar('Limits', [cmin, cmax]); % 固定颜色范围
    
    % subplot(1, 2, 2)
    % surf(x, y, mod(angle(U), 2*pi), 'EdgeColor', 'none', 'FaceAlpha', 0.8)
    % xlabel("x")
    % ylabel("y")
    % zlabel("\phi")
    % title(name+" 相位")
    % cmin = 0; % 数据最小值
    % cmax = 2 * pi; % 数据最大值
    % colorbar('Limits', [cmin, cmax]); % 固定颜色范围
end
