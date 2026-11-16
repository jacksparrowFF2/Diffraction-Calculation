function [pupil,info]=corePlaneAperture(cfg,xx,yy,pupil)
%corePlaneAperture 平面波在物面上的均匀矩形/圆形照明，与几何光瞳取交集。
% 单波长内部空间相干；不同波长在强度上叠加。此处不模拟空间非相干光。
    source=coreSourceDefaults();
    if isfield(cfg,'source')
        for name=string(fieldnames(cfg.source)).', source.(name)=cfg.source.(name); end
    end
    info=struct('enabled',false,'shape',"none",'centerXY',source.centerXY, ...
        'rectangleSize',source.rectangleSize,'circleDiameter',source.circleDiameter, ...
        'requestedArea',NaN,'illuminatedArea',sum(pupil(:).^2)*(cfg.pitch/cfg.samplesPerPeriod)^2, ...
        'coverageFraction',1,'extendsBeyondGrid',false, ...
        'coherence',"每个源模式内部空间相干；不同源模式/波长强度叠加");
    if source.type~="plane", return; end
    mustBeMember(source.apertureShape,["none","rectangle","circle"]);
    if source.apertureShape=="none", return; end
    validateattributes(source.centerXY,{'numeric'},{'size',[1,2],'real','finite'});
    x=xx-source.centerXY(1); y=yy-source.centerXY(2);
    if source.apertureShape=="rectangle"
        validateattributes(source.rectangleSize,{'numeric'},{'size',[1,2],'real','finite','positive'});
        half=source.rectangleSize/2;
        % 半开区间使整像素宽度的矩形不会重复包含两个端点。
        inside=x>=-half(1) & x<half(1) & y>=-half(2) & y<half(2);
        requestedArea=prod(source.rectangleSize);
    else
        validateattributes(source.circleDiameter,{'numeric'},{'scalar','real','finite','positive'});
        half=[1,1]*source.circleDiameter/2;
        inside=x.^2+y.^2<=half(1)^2;
        requestedArea=pi*half(1)^2;
    end
    original=sum(pupil(:).^2);
    pupil=pupil.*inside;
    assert(any(pupil(:)>0),'COEReflection:EmptyIllumination', ...
        '平面波照明与输入物面/光瞳无交集，请检查尺寸和centerXY。');
    dx=cfg.pitch/cfg.samplesPerPeriod;
    lower=[min(xx,[],'all'),min(yy,[],'all')]-dx/2;
    upper=[max(xx,[],'all'),max(yy,[],'all')]+dx/2;
    info.enabled=true; info.shape=source.apertureShape;
    info.requestedArea=requestedArea;
    info.illuminatedArea=sum(pupil(:).^2)*dx^2;
    info.coverageFraction=sum(pupil(:).^2)/original;
    info.extendsBeyondGrid=any(source.centerXY-half<lower-dx*1e-9 ...
        | source.centerXY+half>upper+dx*1e-9);
    if info.extendsBeyondGrid
        fprintf('平面波%s照明超出输入物面：仅计算交集 %.6g mm²，不自动扩大或平铺图样。\n', ...
            info.shape,info.illuminatedArea*1e6);
    end
end
