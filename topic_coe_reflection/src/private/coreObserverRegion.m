function region=coreObserverRegion(cfg,reference)
%coreObserverRegion 镜面零级的源-入瞳光线对贡献区域；仅诊断，不作为衍射硬孔径。
    a=cfg.source.distance; b=cfg.observer.distance;
    rs=b/(a+b)*cfg.source.diameter/2; rp=a/(a+b)*cfg.observer.pupilDiameter/2;
    radius=rs+rp; dx=reference.x(2)-reference.x(1);
    bounds=[reference.x(1)-dx/2,reference.x(end)+dx/2,reference.y(1)-dx/2,reference.y(end)+dx/2];
    extent=1.1*max([radius,abs(bounds)]); coordinate=linspace(-extent,extent,401);
    step=coordinate(2)-coordinate(1); [x,y]=meshgrid(coordinate);
    r=hypot(x,y); density=zeros(size(r));
    if radius==0
        density(201,201)=1/step^2;
    elseif min(rs,rp)==0
        density(r<=radius)=1/(pi*radius^2);
    else
        small=min(rs,rp); large=max(rs,rp);
        overlap=zeros(size(r)); overlap(r<=large-small)=pi*small^2;
        active=r>large-small & r<large+small; d=r(active);
        c1=min(1,max(-1,(d.^2+small^2-large^2)./(2*d*small)));
        c2=min(1,max(-1,(d.^2+large^2-small^2)./(2*d*large)));
        overlap(active)=small^2*acos(c1)+large^2*acos(c2) ...
            -0.5*sqrt(max(0,(-d+small+large).*(d+small-large).*(d-small+large).*(d+small+large)));
        density=overlap/(pi*small^2*pi*large^2);
    end
    illuminated=interp2(reference.x,reference.y,double(reference.pupil>0),x,y,'nearest',0)>0;
    illuminated=illuminated & x>=bounds(1) & x<bounds(2) & y>=bounds(3) & y<bounds(4);
    support=r<=radius;
    if radius==0, support=density>0; end
    inModel=support & illuminated;
    region=struct('coordinate',coordinate,'sourceProjectedDiameter',2*rs, ...
        'pupilProjectedDiameter',2*rp,'diameter',2*radius, ...
        'nominalPeriods',2*radius/cfg.pitch,'geometricArea',pi*radius^2, ...
        'rayPairDensity',density,'support',support,'modeledSupport',inModel, ...
        'modeledArea',sum(inModel,'all')*step^2,'inputBounds',bounds, ...
        'rayWeightInModel',sum(density.*illuminated,'all')/max(sum(density,'all'),realmin), ...
        'rayPairMass',sum(density,'all')*step^2,'diagnosticStep',step, ...
        'meaning',"居中镜面零级；源与入瞳的投影圆盘卷积；标称Pitch数不是相干周期数；未裁切衍射场");
end
