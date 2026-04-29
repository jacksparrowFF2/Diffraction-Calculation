% 输入（匿名函数，物面尺寸xy，物面xy方向分辨率，波长，传播距离，像面尺寸xy，像面xy方向分辨率）
function [Uc1, Ud1, Ud] = fraunhofer_fft_SPH(Uc, xmin, xmax, ymin, ymax, m, n, lambda, z, Xmin, Xmax, Ymin, Ymax, M, N,WaveVector,flag)
    %FRAUNHOFER_FFT 计算连续分布在特定取样条件下的夫琅禾费衍射积分（FFT版）
    %   小写字母为输入光场采样参数，大写字母为输出光场采样参数
    %   lambda 波长
    %   z 传播距离
    %  flag 控制离散化是否开启
    
    % WaveVector = 2 * pi / lambda;
    
    % 如果输入的是匿名函数，则需要进行离散化，flag为1
    if flag == 1
        % 对复合匿名函数进行离散化
        Ud = discretize(Uc, xmin, xmax, ymin, ymax, m, n);
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
    Ud1 = myFFT2(Ud, xmin, xmax, ymin, ymax, Xmin*t, Xmax*t, Ymin*t, Ymax*t, M, N);
    % 此处计算fraunhofer积分中的相位因子,来自于公式
    Ud0 = exp(1j*WaveVector*z+1j*WaveVector/(2 * z)*(X.^2 + Y.^2)) / (1j * lambda * z);
    Ud1 = Ud1 .* Ud0;
    % 将衍射积分计算结果投影到像面
    Uc1 = interpolate(Ud1, Xmin, Xmax, Ymin, Ymax);
end
