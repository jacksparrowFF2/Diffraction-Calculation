%% 验证圆孔的夫琅禾费衍射公式
clear, clc, close all
%% 模型参量
lambdalist = [620e-9,550e-9,450e-9];
% 传播距离,模拟人眼视距2~40cm
z = 0.3;
% 设定周期尺寸,直观单位微米
Pitch = 55e-6;
% 圆孔半径
r = 10e-6;
% 定义物面的光栅类型
% flag=0代表使用相同离散矩阵（自定义的形状或者图片）
% flag=1代表使用不相同的离散矩阵
% flag=2代表使用连续函数
Gratingflag = 1;
% 定义是否仿真像面的衍射效果，平面波为0，图片为1
Lightflag = 1;
% 分辨率，设定单个周期的分辨率
nn = 1000;
% 设定物面X、Y方向周期
Nx_period = 2;         
Ny_period = 2;         
% 设定目标像面单边尺寸，直观单位mm
Fx =40e-3;
Fy =40e-3;
% 设定物体在像面的尺寸，直观单位mm
F_size_x = 20e-3;   % F在像面上的物理宽度
F_size_y = 20e-3;   % F在像面上的物理高度
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
xmin = -Nx_period*Pitch;
% 物面尺寸x方向总最大值
xmax = -xmin;
% 物面尺寸y方向总最小值
ymin = -Ny_period*Pitch;
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

if Lightflag ==1
    % 导入要被衍射的图片
    Input = imread('.\input\USAF500.png');
    U0 = double(rgb2gray(Input));
    figure(1);hold on;
    subplot(1,2,1);
    imshow(U0);
    title('初始振幅'); 
    axis on;axis image;xlabel('x');ylabel('y')
else
    U0 = 1; %不参与最终计算
end


if Gratingflag ==0
    % 创建离散的单元函数
    SinglePattern = CustomPattern(0,Xc,Yc,r,r,0,0);
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
    a1 = r/0.8; b1=r;
    SinglePattern1 = CustomPattern(30,Xc,Yc,a1,b1,0,0);
    % 创建离散的单元函数2，长半轴为a，短半轴为b
    a2 = r; b2=r;
    SinglePattern2 = CustomPattern(0,Xc,Yc,a2,b2,0,0);
    % 创建离散的单元函数3，长半轴为a，短半轴为b
    a3 = r; b3=r;
    SinglePattern3 = CustomPattern(0,Xc,Yc,a3,b3,0,0);
    % 创建离散的单元函数4，长半轴为a，短半轴为b
    a4 = r/0.8; b4=r;
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



% circle2 = @(x, y) (x-10e-6).^2 + (y-10e-6).^2 <= r^2;
% circle = @(x, y)circle1(x,y)&circle2(x,y);

for c = 1:3
    lambda = lambdalist(c);
    % 计算波矢量
    k = 2 * pi / lambda;
    % 数值解-优先计算数值解
    [~, U2] = fraunhofer_fft_SPH(Mask, xmin, xmax, ymin, ymax, Tnn, Tnn, ...
                            lambda, z, Xmin, Xmax, Ymin, Ymax, Tnn, Tnn, ...
                            k,Gratingflag);
    % 绘制图像
    % 计算振幅
    A = abs(U2);
    % 强度分布
    I =A.^2;
    % 对数压缩亮部
    alpha = 300;
    I_log = log(1 + alpha * I);
    I_log = I_log ./ max(I_log(:));
    % 伽马增强暗部
    gamma = 0.25;
    I_enhanced = I_log .^ gamma;

    figure;subplot(1, 2, 1)
    imagesc(X(1,:), Y(:,1), A)
    xlabel("x");ylabel("y");zlabel("|U|")
    axis image;title("光栅振幅图像"+num2str(lambda*1e9) +'nm')
    colormap gray;
    
    subplot(1, 2, 2)
    imagesc(X(1,:), Y(:,1), I_enhanced)
    xlabel("x");ylabel("y");zlabel("|U|")
    axis image;title("光栅强度分布"+num2str(lambda*1e9) +'nm')
    colormap gray;
    
    % RGB(:,:,c) = (I./max(max(I))).^gamma;
    RGB(:,:,c) = I_enhanced;
end

%% RGB合成显示色分离
figure;
imagesc(X(1,:)*1e3, Y(:,1)*1e3, RGB);
axis image;
xlabel('x / mm');
ylabel('y / mm');
title('RGB Far-field Color Separation');
set(gca, 'YDir', 'normal');

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
end 


%% 雾度计算
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

%% 绘图函数
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