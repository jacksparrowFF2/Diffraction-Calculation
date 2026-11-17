function view=coreEyeImage(cfg)
%coreEyeImage 逐非相干源点：入射复场→反射→Fresnel→圆瞳/薄透镜→像面→强度和。
% 反射光路展开为正向传播；法向标量局部反射系数。全局传播相位在取强度前可略。
% 不施加a/(a+b)几何缩放、不卷积源强度或Airy强度、不按图像窗口补能量。
    validateattributes(cfg.periods,{'numeric'},{'size',[1,2],'integer','positive'});
    assert(all(mod(cfg.periods,2)==0),'COEReflection:EyePeriods','图片阵列每轴周期数须为偶数。');
    validateattributes(cfg.samplesPerPeriod,{'numeric'},{'scalar','integer','>=',4});
    validateattributes(cfg.source.distance,{'numeric'},{'scalar','finite','positive'});
    validateattributes(cfg.observer.distance,{'numeric'},{'scalar','finite','positive'});
    validateattributes(cfg.observer.pupilDiameter,{'numeric'},{'scalar','finite','positive'});
    validateattributes(cfg.observer.pupilSamples,{'numeric'},{'scalar','integer','>=',9,'odd'});
    validateattributes(cfg.outputSize,{'numeric'},{'scalar','integer','>=',3,'odd'});
    validateattributes(cfg.observer.angleHalfRangeDeg,{'numeric'},{'scalar','finite','>',0,'<=',10});
    validateattributes(cfg.observer.retinaDistance,{'numeric'},{'scalar','finite','positive'});
    validateattributes(cfg.wavelengths,{'numeric'},{'row','finite','positive','nonempty'});
    validateattributes(cfg.spectralWeights,{'numeric'},{'vector','finite','nonnegative'});
    validateattributes(cfg.channelResponse,{'numeric'},{'2d','finite','nonnegative'});
    validateattributes(cfg.heightMap,{'numeric'},{'real','finite'});
    assert(numel(cfg.spectralWeights)==numel(cfg.wavelengths) && size(cfg.channelResponse,1)==numel(cfg.wavelengths), ...
        'COEReflection:EyeSpectrum','光谱权重/通道维度不符。');
    assert(cfg.mode=="reflection" && cfg.source.type=="point" && cfg.source.integrationMethod=="direct" ...
        && any(cfg.source.wavefront==["paraxial","spherical"]), ...
        'COEReflection:EyeSource','完整成像链要求point/direct及paraxial或spherical入射复场。');
    mustBeMember(cfg.observer.focusMode,["virtualSource","screen","distance"]);
    a=cfg.source.distance; b=cfg.observer.distance;
    focusDistance=a+b;
    if cfg.observer.focusMode=="screen", focusDistance=b; end
    if cfg.observer.focusMode=="distance", focusDistance=cfg.observer.focusDistance; end
    validateattributes(focusDistance,{'numeric'},{'scalar','finite','positive'});
    dx=cfg.pitch/cfg.samplesPerPeriod;
    shape=cfg.periods([2,1])*cfg.samplesPerPeriod;
    x=((0:shape(2)-1)-floor(shape(2)/2))*dx;
    y=((0:shape(1)-1)-floor(shape(1)/2))*dx;
    [xx,yy]=meshgrid(x,y); aperture=double(hypot(xx,yy)<=cfg.pupilRadius);
    assert(any(aperture,'all'),'COEReflection:EyeAperture','计算面上无受照区域。');
    p=linspace(-cfg.observer.pupilDiameter/2,cfg.observer.pupilDiameter/2,cfg.observer.pupilSamples);
    dp=p(2)-p(1); [px,py]=meshgrid(p); pupil=double(hypot(px,py)<=cfg.observer.pupilDiameter/2);
    theta=linspace(-deg2rad(cfg.observer.angleHalfRangeDeg),deg2rad(cfg.observer.angleHalfRangeDeg),cfg.outputSize);
    dt=theta(2)-theta(1); wavelengths=cfg.wavelengths; nLambda=numel(wavelengths);
    % 除输出角域外，检查瞳面二次相位和物面有限支撑对应的场带宽。
    pupilFrequency=max(abs(theta))/min(wavelengths);
    fieldFrequency=max(abs([x,y]))/(min(wavelengths)*b) ...
        +max(abs(p))*abs(1/b-1/focusDistance)/min(wavelengths);
    objectFrequency=(max(abs([x,y]))+norm(cfg.source.centerXY)+cfg.source.diameter/2)/(min(wavelengths)*a) ...
        +(max(abs([x,y]))+max(abs(p)))/(min(wavelengths)*b);
    assert(objectFrequency<1/(2*dx),'COEReflection:EyeObjectSampling', ...
        '物面复场采样不足：提高samplesPerPeriod或缩小物面/源/瞳范围。');
    assert(max(pupilFrequency,fieldFrequency)<1/(2*dp), ...
        'COEReflection:EyePupilSampling','瞳面采样不足：提高observer.pupilSamples或缩小观察角域。');
    assert(dt<=min(wavelengths)/(2*cfg.observer.pupilDiameter), ...
        'COEReflection:EyeImageSampling','像面强度采样不足：提高outputSize或缩小观察角域。');
    incidentArea=sum(aperture,'all')*dx^2;
    energy=zeros(cfg.outputSize,cfg.outputSize,nLambda);
    baselineEnergy=energy; pupilPower=zeros(1,nLambda); reflectedPower=pupilPower;
    baselinePupilPower=pupilPower; pupilIntensity=zeros(size(pupil));
    hasHeight=any(cfg.heightMap(:)~=0); stateCount=1+double(hasHeight);
    started=tic;
    for wave=1:nLambda
        lambda=wavelengths(wave);
        [incoming,weights,sourceInfo]=coreBuildIllumination(cfg,xx,yy,aperture,lambda);
        coefficient=reflectionCoefficient(cfg,shape,wave);
        if hasHeight
            assert(isscalar(cfg.heightMap)||isequal(size(cfg.heightMap),shape),'COEReflection:EyeHeight','heightMap尺寸不符。');
            phase=exp(4i*pi*cfg.deltaIndex*cfg.heightMap/lambda);
        else
            phase=1;
        end
        objectChirp=exp(1i*pi*(xx.^2+yy.^2)/(lambda*b));
        pupilChirp=exp(1i*pi*(px.^2+py.^2)/(lambda*b));
        % 透镜相位与像距传播二次相位合并：1/f=1/v+1/focusDistance。
        lensAndRetina=exp(-1i*pi*(px.^2+py.^2)/(lambda*focusDistance));
        for mode=1:numel(weights)
            for state=1:stateCount
                object=incoming{mode}.*aperture.*coefficient;
                if state==stateCount, object=object.*phase; end
                atPupil=coreFourierIntegral(object.*objectChirp,x,y,p/(lambda*b),p/(lambda*b)) ...
                    .*pupilChirp/(1i*lambda*b);
                admitted=atPupil.*pupil;
                image=coreFourierIntegral(admitted.*lensAndRetina,p,p,theta/lambda,theta/lambda)/lambda;
                contribution=weights(mode)*abs(image).^2*dt^2/incidentArea;
                collected=weights(mode)*sum(abs(admitted).^2,'all')*dp^2/incidentArea;
                if state==stateCount
                    energy(:,:,wave)=energy(:,:,wave)+contribution;
                    pupilPower(wave)=pupilPower(wave)+collected;
                    reflectedPower(wave)=reflectedPower(wave)+weights(mode)*sum(abs(object).^2,'all')*dx^2/incidentArea;
                    if wave==ceil(nLambda/2), pupilIntensity=pupilIntensity+weights(mode)*abs(admitted).^2; end
                else
                    baselineEnergy(:,:,wave)=baselineEnergy(:,:,wave)+contribution;
                    baselinePupilPower(wave)=baselinePupilPower(wave)+collected;
                end
            end
            if mode==1 || mod(mode,20)==0 || mode==numel(weights)
                fprintf('复场成像 %.0f nm：源点 %d/%d，已用%.1f s\n',lambda*1e9,mode,numel(weights),toc(started));
                drawnow limitrate;
            end
        end
    end
    if ~hasHeight, baselineEnergy=energy; baselinePupilPower=pupilPower; end
    channelWeights=cfg.channelResponse.*cfg.spectralWeights(:);
    view=struct('model',"Per-source coherent Fresnel pupil imaging; incoherent source intensity sum", ...
        'config',cfg,'coordinate',theta,'coordinateUnit',"paraxial image angle (retina position / retina distance), radian", ...
        'retinaCoordinate',theta*cfg.observer.retinaDistance,'spectralEnergy',energy, ...
        'baselineSpectralEnergy',baselineEnergy,'binSolidAngle',dt^2,'heightPhaseEnabled',hasHeight, ...
        'pupilPowerFraction',pupilPower,'baselinePupilPowerFraction',baselinePupilPower, ...
        'reflectedPowerFraction',reflectedPower,'focusDistance',focusDistance, ...
        'lensFocalLength',1/(1/cfg.observer.retinaDistance+1/focusDistance), ...
        'pupilCoordinate',p,'pupilIrradianceRatio',pupilIntensity,'sourceInfo',sourceInfo, ...
        'sourceWeights',weights,'incidentArea',incidentArea,'elapsedSeconds',toc(started));
    view.channelSignal=reshape(reshape(energy,[],nLambda)*channelWeights,cfg.outputSize,cfg.outputSize,[]);
    view.baselineChannelSignal=reshape(reshape(baselineEnergy,[],nLambda)*channelWeights,cfg.outputSize,cfg.outputSize,[]);
    view.windowPowerFraction=reshape(sum(sum(energy,1),2),1,[]);
    view.windowToPupilFraction=view.windowPowerFraction./max(pupilPower,realmin);
    view.assumptions="标量近轴Fresnel；薄透镜无像差；居中圆瞳；独立源点；同波长单源保留全孔径相干；像素中心密度求积";
end

function coefficient=reflectionCoefficient(cfg,shape,wave)
    if ~isempty(cfg.complexCoefficient)
        assert(size(cfg.complexCoefficient,3)==1||size(cfg.complexCoefficient,3)==numel(cfg.wavelengths), ...
            'COEReflection:EyeCoefficient','复反射系数波长维度不符。');
        assert(isscalar(cfg.complexCoefficient)||isequal(size(cfg.complexCoefficient,1:2),shape), ...
            'COEReflection:EyeCoefficient','复反射系数尺寸不符。');
        coefficient=cfg.complexCoefficient(:,:,min(wave,size(cfg.complexCoefficient,3))).*ones(shape);
    else
        assert(isfield(cfg,'coeStack')&&cfg.coeStack.enabled, ...
            'COEReflection:EyeCoefficient','完整成像链要求图片膜堆或显式complexCoefficient。');
        values=cfg.coeStack.coefficientLookup(:,wave);
        coefficient=reshape(values(double(cfg.coeStack.regionMap(:))+1),shape);
    end
    assert(all(isfinite(coefficient),'all')&&all(abs(coefficient)<=1+1e-12,'all'), ...
        'COEReflection:EyeCoefficient','复反射系数应有限且满足被动反射。');
end
