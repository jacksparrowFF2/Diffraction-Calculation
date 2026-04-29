%% 验证圆孔的夫琅禾费衍射公式
clear, clc, close all

%% 模型参量
lambda = 500e-9;
% 传播距离,模拟人眼视距2~40cm
z = 0.3;
% 圆孔半径
r = 10e-6;
% 初始振幅
U0 = 1;
% 定义物面类型。flag
flag = 0;
% 分辨率，设定单个周期的分辨率
nn = 100;
% 设定物面X、Y方向周期
Nx_period = 3;         
Ny_period = 3;         
% 设定周期尺寸
Pitch = 55e-6;
% 设定目标像面尺寸，直观单位mm
Fx =20e-3;
Fy =20e-3;
%设定目标角，用于寻找指定角度之外和之内的强度和，单位是°，黄度是2.5°
Angle = 2.5;
% 计算目标角对应传播距离下的长度
FocusR = tand(Angle)*z;
fprintf('像面上2.5°范围半径为：%g mm\n', FocusR*1e3);


% 导入光源信息（光源即为要被评估经过衍射后的图案）
RGB = imread('.\input\输入图片2.jpg');
InputPic = double(rgb2gray(RGB));

% 设定单个周期的最大与最小值
x_cell_min = -Pitch/2;
x_cell_max = Pitch/2;
y_cell_min = -Pitch/2;
y_cell_max = Pitch/2;
% 设定单个周期的格点坐标
x_cell = linspace(x_cell_min,x_cell_max,nn);
y_cell = linspace(y_cell_min,y_cell_max,nn);

[Xc, Yc] = meshgrid(x_cell, y_cell);





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

% 计算波矢量
k = 2 * pi / lambda;
if flag ==1
    % 创建用于作为物面的函数，一定是连续函数
    Mask = @(Xc, Yc) Xc.^2 + Yc.^2 <= r^2;
else
    % 创建离散的单元函数
    SinglePattern = double(Xc.^2 + Yc.^2 <= r^2);
    % 在x和y方向上分别重复Nx_period、Ny_period次周期
    Mask = repmat(SinglePattern, Ny_period, Nx_period);
    figure;
    imshow(Mask, []);
    title('单个Pitch内形状: Circular Aperture'); 
    axis on;axis image;xlabel('x');ylabel('y')
end



% circle2 = @(x, y) (x-10e-6).^2 + (y-10e-6).^2 <= r^2;
% circle = @(x, y)circle1(x,y)&circle2(x,y);

% 创建用于作为屋面的离散Map
%% 数值解-优先计算数值解
% [~, U2] = fraunhofer_fft(@(x, y)x.^2+y.^2 <= r^2, xmin, xmax, ymin, ymax, nn, nn, lambda, z, Xmin, Xmax, Ymin, Ymax, nn, nn);
[~, U2] = fraunhofer_fft_SPH(Mask, xmin, xmax, ymin, ymax, Tnn, Tnn, ...
                        lambda, z, Xmin, Xmax, Ymin, Ymax, Tnn, Tnn, ...
                        k,flag);
work(X, Y, U2, "数值解");

% 计算振幅
A = abs(U2);
% 强度分布
I =A.^2;

% 构建像面上在指定角度下的范围遮罩
AngleMask = double(X.^2 + Y.^2 <= FocusR^2);
figure;
% imshow(AngleMask, []);
imagesc(X(1,:),Y(:,1),AngleMask)
title('像面角度遮罩: Circular Aperture'); 
axis on;axis image;xlabel('x');ylabel('y')

%计算总能量（因为是计算比例，所以没有乘以像面xy步长）
I_Total = sum(I,'all');
fprintf('总能量：%g\n', I_Total);
% 计算2.5度范围内的能量
I_Angle = sum(I.*AngleMask,'all');
fprintf('总能量：%g\n', I_Angle);
%计算雾度
WD = 1-I_Angle/I_Total;
fprintf('理论雾度：%g %%\n', WD*100);

%% 理论解-圆孔衍射的理论解是一阶贝塞尔函数
% U1 = @(x, y)exp(1j*k*(z + (x.^2 + y.^2) / (2 * z))) * r ./ (1j * sqrt(x.^2+y.^2)) .* besselj(1, 2*pi*r*sqrt(x.^2+y.^2)/(lambda * z));
% U1 = U1(X, Y);
% work(X, Y, U1, "理论解");
%% 误差函数
% work(X, Y, U2-U1, "数值-理论");

%% 绘图函数
function work(x, y, U, name)

% 振幅
A = abs(U);
% 强度
I =A.^2;
figure
subplot(1, 2, 1)
% surf(x, y, abs(U), 'EdgeColor', 'none', 'FaceAlpha', 0.8)
imagesc(x(1,:), y(:,1), A)
xlabel("x")
ylabel("y")
zlabel("|U|")
axis equal
title(name+" 振幅图像")
cmin = min(A(:)); % 数据最小值
cmax = max(A(:)); % 数据最大值
colorbar('Limits', [cmin, cmax]); % 固定颜色范围


subplot(1, 2, 2)
imagesc(x(1,:), y(:,1), I)
xlabel("x")
ylabel("y")
zlabel("|U|^2")
axis equal
title(name+" 强度分布")
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