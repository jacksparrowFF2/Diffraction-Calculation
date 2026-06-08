function Csame = fftconv2_same(A, B)
% 计算二维卷积，返回与 A 同尺寸结果
% B 默认看作中心在矩阵中心的核

    [Ma, Na] = size(A);
    [Mb, Nb] = size(B);

    M = Ma + Mb - 1;
    N = Na + Nb - 1;

    % 卷积核先 ifftshift，使“中心”变成 FFT 所需形式
    % B2 = ifftshift(B);

    FA = fft2(A, M, N);
    FB = fft2(B, M, N);

    C = ifft2(FA .* FB);

    % 截取 same 区域
    rowStart = floor(Mb/2) + 1;
    colStart = floor(Nb/2) + 1;

    Csame = C(rowStart:rowStart+Ma-1, colStart:colStart+Na-1);
end


% function C = fftconv2_same(A, B)
% 
%     sizeA = size(A);
%     sizeB = size(B);
% 
%     outSize = sizeA + sizeB - 1;
% 
%     FA = fft2(A, outSize(1), outSize(2));
%     FB = fft2(B, outSize(1), outSize(2));
% 
%     Cfull = real(ifft2(FA .* FB));
% 
%     start1 = floor(sizeB(1)/2) + 1;
%     start2 = floor(sizeB(2)/2) + 1;
% 
%     C = Cfull(start1:start1+sizeA(1)-1, ...
%               start2:start2+sizeA(2)-1);
% end