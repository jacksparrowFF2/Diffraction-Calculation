function output=coreFourierIntegral(input,x,y,fx,fy)
%coreFourierIntegral 二维复傅里叶积分；保留相位与dx*dy，指定均匀频率网格。
% Bluestein卷积；不要求输出频率落在FFT原生频点，不依赖Signal Processing Toolbox。
    output=columns(columns(input.',x,fx).',y,fy);
end

function output=columns(input,coordinate,frequency)
    n=size(input,1); m=numel(frequency);
    assert(numel(coordinate)==n && n>=2 && m>=2,'COEReflection:FourierShape','傅里叶积分网格尺寸不符。');
    dx=coordinate(2)-coordinate(1); df=frequency(2)-frequency(1);
    assert(dx>0 && df>0,'COEReflection:FourierGrid','傅里叶积分要求递增均匀网格。');
    assert(max(abs(diff(coordinate)-dx))<1e-10*dx && max(abs(diff(frequency)-df))<1e-10*df, ...
        'COEReflection:FourierGrid','傅里叶积分要求均匀网格。');
    ni=(0:n-1).'; mi=(0:m-1).'; count=2^nextpow2(n+m-1);
    kernel=complex(zeros(count,1));
    kernel(1:m)=exp(1i*pi*dx*df*mi.^2);
    kernel(count-n+2:count)=exp(1i*pi*dx*df*((-n+1:-1).').^2);
    chirped=input.*exp(-1i*pi*(dx*df*ni.^2+2*dx*frequency(1)*ni));
    value=ifft(fft(chirped,count,1).*fft(kernel),[],1);
    output=value(1:m,:).*exp(-1i*pi*(dx*df*mi.^2 ...
        +2*df*coordinate(1)*mi+2*coordinate(1)*frequency(1)))*dx;
end
