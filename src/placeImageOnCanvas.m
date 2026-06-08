function canvas = placeImageOnCanvas(img, Nx, Ny, H_obj, W_obj, Fx, Fy)
% img    : 原始字母图像（0~1）
% Ny,Nx  : 仿真画布大小
% H_obj  : 字母图像物理高度 [m]
% W_obj  : 字母图像物理宽度 [m]
% dy,dx  : 仿真网格采样步长 [m]
    

    % 计算仿真网格采样步长[m]
    dx = Fx/Nx;
    dy = Fy/Ny;
    % 计算该物理尺寸在仿真网格中对应多少像素
    Nh = max(1, round(H_obj / dy));
    Nw = max(1, round(W_obj / dx));

    % 重采样到指定步长对应的分辨率
    img_rs = imresize(img, [Nh, Nw], 'bilinear');

    % 创建与像面同分辨率的画布
    canvas = zeros(Ny, Nx);
    % 计算重采样后的图片在画布中的起始点
    y1 = floor((Ny - Nh)/2) + 1;
    x1 = floor((Nx - Nw)/2) + 1;
    % 计算重采样后的图片在画布中的结束点
    y2 = y1 + Nh - 1;
    x2 = x1 + Nw - 1;
    % 若超出边界则截断，获得校正后的起始和结束点
    yy1 = max(1, y1); yy2 = min(Ny, y2);
    xx1 = max(1, x1); xx2 = min(Nx, x2);
    % 计算重采样后的图片能被放入画布的区域
    sy1 = 1 + (yy1 - y1);
    sy2 = Nh - (y2 - yy2);
    sx1 = 1 + (xx1 - x1);
    sx2 = Nw - (x2 - xx2);
    % 在画布中赋予
    canvas(yy1:yy2, xx1:xx2) = img_rs(sy1:sy2, sx1:sx2);
end