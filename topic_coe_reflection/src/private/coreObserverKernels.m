function detail=coreObserverKernels(cfg,step)
%coreObserverKernels 非相干圆盘源角核及理想圆形入瞳强度PSF。
% 角图强度卷积是额外成像近似，不把瞳孔内的点当作独立光源。
    sourceRadius=cfg.source.diameter/(2*(cfg.source.distance+cfg.observer.distance));
    sourceHalf=ceil(sourceRadius/step);
    sourceAxis=(-sourceHalf:sourceHalf)*step; [x,y]=meshgrid(sourceAxis);
    sourceKernel=zeros(size(x));
    if sourceRadius==0
        sourceKernel=1;
    else
        for u=(-2:2)*step/5
            for v=(-2:2)*step/5
                sourceKernel=sourceKernel+double(hypot(x+u,y+v)<=sourceRadius);
            end
        end
        sourceKernel=sourceKernel/sum(sourceKernel,'all'); % 源点概率积分为1。
    end
    wavelengths=cfg.wavelengths;
    pupilKernels=cell(1,numel(wavelengths)); retained=ones(size(wavelengths));
    firstZero=NaN(size(wavelengths)); pupilHalf=0;
    for k=1:numel(wavelengths)
        pupilKernels{k}=1;
        if cfg.observer.pupilModel~="airyIntensity" || cfg.observer.pupilDiameter==0, continue; end
        firstZero(k)=1.2196698912665045*wavelengths(k)/cfg.observer.pupilDiameter;
        radius=cfg.observer.airyRadiusFactor*firstZero(k);
        half=ceil(radius/step); pupilHalf=max(pupilHalf,half);
        axis=(-half:half)*step; [px,py]=meshgrid(axis); kernel=zeros(size(px));
        for u=(-1:1)*step/3
            for v=(-1:1)*step/3
                r=hypot(px+u,py+v); z=pi*cfg.observer.pupilDiameter*r/wavelengths(k);
                amplitude=ones(size(z)); active=z~=0;
                amplitude(active)=2*besselj(1,z(active))./z(active);
                kernel=kernel+amplitude.^2.*(r<=radius);
            end
        end
        edge=pi*cfg.observer.pupilDiameter*radius/wavelengths(k);
        retained(k)=1-besselj(0,edge)^2-besselj(1,edge)^2;
        % 校准到解析包围能量，不把截尾核重标定为1。
        pupilKernels{k}=kernel*(retained(k)/sum(kernel,'all'));
    end
    detail=struct('sourceKernel',sourceKernel,'pupilKernels',{pupilKernels}, ...
        'stepRad',step,'paddingBins',sourceHalf+pupilHalf, ...
        'sourceAngularDiameterRad',2*sourceRadius,'airyFirstZeroRad',firstZero, ...
        'pupilRetainedMass',retained,'sourceKernelMass',sum(sourceKernel,'all'), ...
        'pupilModel',cfg.observer.pupilModel, ...
        'assumptions',"均匀非相干圆盘源；局部移不变角图；airyIntensity额外采用非相干强度PSF，无像差/离焦/相干像面传播");
end
