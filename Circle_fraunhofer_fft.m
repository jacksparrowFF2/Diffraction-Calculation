%% 验证圆孔的夫琅禾费衍射公式
clear, clc, close all

%% 模型参量
lambda = 500e-9;
% 传播距离,模拟人眼视距2~40cm
z = 0.3;
% 圆孔半径
r = 20e-6;
% 初始振幅
U0 = 1;
% 分辨率，设定单个周期的分辨率
nn = 100;
% 设定物面X、Y方向周期
Nx_period = 1;         % x方向周期数
Ny_period = 1;         % y方向周期数
% 设定周期尺寸
Pitch = 50e-6;
% 设定单个周期的最大与最小值
x_cell_min = -Pitch/2;
x_cell_max = Pitch/2;
y_cell_min = -Pitch/2;
y_cell_max = Pitch/2;
x_cell = linspace(x_cell_min,x_cell_max,nn);
y_cell = linspace(y_cell_min,y_cell_max,nn);

[Xc, Yc] = meshgrid(x_cell, y_cell);

cellPattern = double(Xc.^2 + Yc.^2 <= r^2);

figure;
imshow(cellPattern, []);
title('Single Unit Cell: Circular Aperture'); 
axis on;axis image;xlabel('x');ylabel('y')


% 像面与物面尺寸之比
ImageMaskRatio = 1e3;
% 物面尺寸x方向最小值
xmin = -2*r;
% 物面尺寸x方向最大值
xmax = -xmin;
% 物面尺寸y方向最小值
ymin = -2*r;
% 物面尺寸y方向最大值
ymax = -ymin;
% 像面尺寸x方向最小值
Xmin = xmin*ImageMaskRatio;
% 像面尺寸x方向最大值
Xmax = -Xmin;
% 像面尺寸y方向最小值
Ymin = ymin*ImageMaskRatio;
% 像面尺寸y方向最大值
Ymax = -Ymin;

% 输出像面尺寸
fprintf('像面尺寸为：%g x %g mm\n', Xmin*2*1e3, Ymin*2*1e3);
% 构建物面x方向矢量
x = linspace(xmin, xmax, nn);
% 构建物面y方向矢量
y = linspace(ymin, ymax, nn);
% 构建物面xy离散坐标点
[x, y] = meshgrid(x, y);
% 构建像面X方向矢量
X = linspace(Xmin, Xmax, nn);
% 构建像面Y方向矢量
Y = linspace(Ymin, Ymax, nn);
% 构建像面XY离散坐标点
[X, Y] = meshgrid(X, Y);

% 计算波矢量
k = 2 * pi / lambda;
% 创建函数
circle = @(x, y) x.^2 + y.^2 <= r^2;

% circle2 = @(x, y) (x-10e-6).^2 + (y-10e-6).^2 <= r^2;
% circle = @(x, y)circle1(x,y)&circle2(x,y);
%% 数值解-优先计算数值解
% [~, U2] = fraunhofer_fft(@(x, y)x.^2+y.^2 <= r^2, xmin, xmax, ymin, ymax, nn, nn, lambda, z, Xmin, Xmax, Ymin, Ymax, nn, nn);
[~, U2] = fraunhofer_fft_SPH(circle, xmin, xmax, ymin, ymax, nn, nn, ...
                        lambda, z, Xmin, Xmax, Ymin, Ymax, nn, nn, ...
                        k,1);
work(X, Y, U2, "数值解");


%% 理论解-圆孔衍射的理论解是一阶贝塞尔函数
% U1 = @(x, y)exp(1j*k*(z + (x.^2 + y.^2) / (2 * z))) * r ./ (1j * sqrt(x.^2+y.^2)) .* besselj(1, 2*pi*r*sqrt(x.^2+y.^2)/(lambda * z));
% U1 = U1(X, Y);
% work(X, Y, U1, "理论解");
%% 误差函数
% work(X, Y, U2-U1, "数值-理论");

%% 绘图函数
function work(x, y, U, name)
figure
% subplot(1, 2, 1)
% surf(x, y, abs(U), 'EdgeColor', 'none', 'FaceAlpha', 0.8)
imagesc(x(1,:), y(:,1), abs(U))
xlabel("x")
ylabel("y")
zlabel("|U|")
% axis equal
title(name+" 振幅")
cmin = min(abs(U(:))); % 数据最小值
cmax = max(abs(U(:))); % 数据最大值
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