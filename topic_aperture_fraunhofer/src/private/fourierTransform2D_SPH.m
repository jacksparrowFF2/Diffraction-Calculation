function Ud1 = fourierTransform2D_SPH(Ud, xmin, xmax, ymin, ymax, Xmin, Xmax, Ymin, Ymax, M, N)
%FOURIERTRANSFORM2D_SPH 在指定空间频率网格上计算二维傅里叶变换。
%   _SPH 表示“曾经理解过并添加了中文注释”。
%   作用：对输入复振幅 Ud 做二维傅里叶积分的离散近似；不包含光学传播。
%   输入 xmin/xmax/ymin/ymax 是物面坐标范围；大写 Xmin/Xmax/Ymin/Ymax
%   在本函数中是空间频率范围（物面单位为 m 时，频率单位为 1/m）。
%   M 为输出列数，N 为输出行数；输出 Ud1 是 N×M 的复数变换结果。
%   与 fft2 的区别：输出频率范围和采样点数可自行指定，但仍须检查采样精度。
%   方法：用二次相位因子将傅里叶求和改写为卷积，再用 FFT 加速；
%   最后乘 kx*ky，计入物面积分面积元。
%   调用关系：主程序 -> fraunhoferDiffraction_SPH -> 本函数。
%   分工：本函数只算傅里叶积分；上层函数负责波长、距离及传播相位。
    %MYFFT2 求离散分布Ud(X0,Y0)对应的连续分布Uc(x0,y0)的傅里叶变换Uc1(x,y)=∫∫Uc(x0,y0)exp(-i2pi(x0*x+y0*y))dx0dy0对应的离散分布Ud1(X,Y)
    %   xmin,xmax,ymin,ymax 原函数的积分范围
    %   Xmin,Xmax,Ymin,Ymax 像函数的计算范围
    %   M,N 像函数的网格数量

    % 获得物面分布的列数
    m = size(Ud, 2);
    % 获得物面分布的行数
    n = size(Ud, 1);
    % 计算物面x方向步长
    kx = (xmax - xmin) / (m - 1);
    % 计算物面y方向步长
    ky = (ymax - ymin) / (n - 1);
    % 计算像面X方向步长
    kX = (Xmax - Xmin) / (M - 1);
    % 计算像面Y方向步长
    kY = (Ymax - Ymin) / (N - 1);
    % 生成物面x坐标序列号
    ix = 0:m - 1;
    % 生成物面y坐标序列号
    iy = 0:n - 1;
    % 生成像面x坐标序列号
    ox = 0:M - 1;
    % 生成像面y坐标序列号
    oy = 0:N - 1;
    % 为卷积构建所有可能的索引差
    IX = [ox, -m:-1];
    IY = [oy, -n:-1];
    % 构建物面索引坐标序列号
    [ix, iy] = meshgrid(ix, iy);
    % 构建像面索引坐标序列号
    [ox, oy] = meshgrid(ox, oy);
    % 构建适用于FFT的循环卷积索引
    [IX, IY] = meshgrid(IX, IY);
    % 计算卷积的第一部分
    U1 = Ud .* exp(-1j*pi*(kx * kX * ix.^2 + ky * kY * iy.^2 + 2 * kx * Xmin * ix + 2 * ky * Ymin * iy));
    % 拓展矩阵
    U1(n+N, m+M) = 0;
    % 计算卷积的第二部分
    U2 = exp(1j*pi*(kx * kX * IX.^2 + ky * kY * IY.^2));
    % 利用卷积定理两个函数的卷积等于两个函数傅里叶变换的逆傅里叶变换
    Ud1 = ifft2(fft2(U1).*fft2(U2)); % 卷积定理
    % Ud1 = cconv2(U1, U2);
    % 乘以相位因子，获得最终矩阵
    Ud1 = Ud1(1:N, 1:M) .* exp(-1j*pi*(kx * kX * ox.^2 + ky * kY * oy.^2 + 2 * kX * xmin * ox + 2 * kY * ymin * oy + 2 * xmin * Xmin + 2 * ymin * Ymin)) * kx * ky;
end
