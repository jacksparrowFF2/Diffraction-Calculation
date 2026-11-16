function [energy,irradiance,clearEnergy,samples,detail] = corePropagateSourceConvolution(field,incoming,x,y,coordinate,lambda,cfg,incidentArea)
%corePropagateSourceConvolution 近轴均匀非相干圆盘源；逐波长强度线性卷积。
% plane：倾斜平面波的平移模型；paraxial：保留中心球面波二次相位。
% 两者均假设各源点物面振幅相同，不包含精确球面距离/辐射包络变化。
    step=coordinate(2)-coordinate(1);
    radius=cfg.source.diameter/(2*cfg.source.distance); % 小角度坐标，rad
    half=ceil(radius/step);
    offsets=(-half:half)*step;
    [sx,sy]=meshgrid(offsets);
    kernel=double(hypot(sx,sy)<=radius+eps(max(radius,step)));
    kernel=kernel/sum(kernel,'all');
    padding=half+1;
    extended=coordinate(1)+(-padding:numel(coordinate)+padding-1)*step;
    assert(max(abs(extended))<=deg2rad(10), ...
        'COEReflection:ConvolutionExtent','卷积补边后的角域超过近轴模型10°范围。');
    [point,~,samples]=corePropagateReceiver(field,x,y,extended,lambda,cfg,incidentArea);
    clearPoint=corePropagateReceiver(incoming,x,y,extended,lambda,cfg,incidentArea);
    [convolved,massError]=linearConvolve(point,kernel);
    clearConvolved=linearConvolve(clearPoint,kernel);
    selected=padding+(1:numel(coordinate));
    energy=convolved(selected,selected);
    clearEnergy=clearConvolved(selected,selected);
    irradiance=energy*incidentArea/(cfg.propagationDistance*step)^2;
    detail=struct('kernel',kernel,'kernelStepRad',step,'kernelMass',sum(kernel,'all'), ...
        'angularRadiusDeg',rad2deg(radius),'paddingBins',padding, ...
        'paddedHalfAngleDeg',rad2deg(max(abs(extended))), ...
        'pointEnergy',point(selected,selected),'pointClearEnergy',clearPoint(selected,selected), ...
        'massError',massError,'wavefront',cfg.source.wavefront, ...
        'assumption',"近轴、相同振幅、均匀非相干圆盘；卷积在强度上执行，保留窗口外向内贡献，不归一化最终裁切结果");
end

function [same,massError]=linearConvolve(value,kernel)
    shape=size(value)+size(kernel)-1;
    fftShape=2.^nextpow2(shape);
    full=real(ifft2(fft2(value,fftShape(1),fftShape(2)).*fft2(kernel,fftShape(1),fftShape(2))));
    full=full(1:shape(1),1:shape(2));
    massError=abs(sum(full,'all')-sum(value,'all'));
    assert(massError<=1e-11*max(1,sum(value,'all')), ...
        'COEReflection:ConvolutionEnergy','线性卷积能量检查失败。');
    full=max(full,0); % FFT舍入级负值；不重标定总功率
    first=(size(kernel)-1)/2;
    same=full(first(1)+(1:size(value,1)),first(2)+(1:size(value,2)));
end
