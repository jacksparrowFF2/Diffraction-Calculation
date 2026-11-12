function [Uc1, Ud1, Ud] = fraunhoferDiffraction_SPH(Uc, xmin, xmax, ymin, ymax, m, n, lambda, z, Xmin, Xmax, Ymin, Ymax, M, N,WaveVector,flag)
%FRAUNHOFERDIFFRACTION_SPH 计算指定观察面网格上的夫琅禾费衍射复振幅。
%   _SPH 表示“曾经理解过并添加了中文注释”。
%   作用：在夫琅禾费远场近似下，把物面复振幅传播到距离 z 的观察面。
%   输入 lambda 为波长，z 为传播距离，WaveVector 应为 2*pi/lambda；
%   大写 Xmin/Xmax/Ymin/Ymax 在本函数中是观察面实际坐标范围，
%   与物面坐标、lambda、z 使用一致的长度单位；M/N 为输出列数/行数。
%   flag==2 时先采样匿名函数 Uc，其它情况直接将 Uc 作为离散矩阵。
%   步骤：将观察面坐标换算为频率 fx=X/(lambda*z)、fy=Y/(lambda*z)，
%   调用 fourierTransform2D_SPH，再乘传播系数及观察面二次相位：
%   U = exp(i*k*z)/(i*lambda*z) * exp(i*k*(X^2+Y^2)/(2*z)) * F(fx,fy)。
%   输出 Uc1 为观察面复振幅的插值函数，Ud1 为 N×M 复振幅矩阵，
%   Ud 为物面离散复振幅；衍射强度需另算 abs(Ud1).^2。
%   分工：本函数负责光学传播；fourierTransform2D_SPH 只负责傅里叶积分。
    %FRAUNHOFER_FFT 计算连续分布在特定取样条件下的夫琅禾费衍射积分（FFT版）
    %   小写字母为输入光场采样参数，大写字母为输出光场采样参数
    %   lambda 波长
    %   z 传播距离
    %  flag 控制离散化是否开启

    % WaveVector = 2 * pi / lambda;

    % 如果输入的是匿名函数，则需要进行离散化，flag为1
    if flag == 2
        % 对复合匿名函数进行离散化
        Ud = Selfdiscretize(Uc, xmin, xmax, ymin, ymax, m, n);
    % 如果输入的是图片或离散化后的Mask,则flag为0,Ud=Uc
    else
        Ud = Uc;
    end

    % 构建像面X方向矢量
    X = linspace(Xmin, Xmax, M);
    % 构建像面Y方向矢量
    Y = linspace(Ymin, Ymax, N);
    % 构建像面XY离散坐标点
    [X, Y] = meshgrid(X, Y);
    % 计算缩放因子，这个缩放因子来自于fraunhofer衍射积分中的相位因子，这样就可以将
    % fraunhofer衍射积分简化成傅里叶变换形式
    t = 1 / (lambda * z);
    Ud1 = fourierTransform2D_SPH(Ud, xmin, xmax, ymin, ymax, Xmin*t, Xmax*t, Ymin*t, Ymax*t, M, N);
    % 此处计算fraunhofer积分中的相位因子,来自于公式
    Ud0 = exp(1j*WaveVector*z+1j*WaveVector/(2 * z)*(X.^2 + Y.^2)) / (1j * lambda * z);
    Ud1 = Ud1 .* Ud0;
    % 将衍射积分计算结果投影到像面
    Uc1 = interpolate(Ud1, Xmin, Xmax, Ymin, Ymax);
end
