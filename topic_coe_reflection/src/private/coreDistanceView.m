function view = coreDistanceView(cfg)
%coreDistanceView theta=a/(a+b)*q；可选有限圆盘源及理想入瞳强度PSF。
% 在原始逐波长能量上积分，接收几何区域独立输出，不裁切反射图样。
    validateattributes(cfg.source.distance,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(cfg.source.diameter,{'numeric'},{'scalar','real','finite','nonnegative'});
    validateattributes(cfg.observer.distance,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(cfg.observer.pupilDiameter,{'numeric'},{'scalar','real','finite','nonnegative'});
    validateattributes(cfg.observer.airyRadiusFactor,{'numeric'},{'scalar','real','finite','>=',1});
    validateattributes(cfg.observer.angleHalfRangeDeg,{'numeric'},{'scalar','real','finite','positive'});
    mustBeMember(cfg.observer.pupilModel,["geometryOnly","airyIntensity"]);
    mustBeMember(cfg.observer.referenceMethod,["quadrature","fftBins"]);
    validateattributes(cfg.observer.referencePaddingFactor,{'numeric'},{'scalar','integer','positive','finite'});
    assert(cfg.mode=="reflection" && cfg.source.type=="point" ...
        && all(cfg.source.centerXY==0) ...
        && cfg.source.integrationMethod=="direct", ...
        'COEReflection:DistanceViewSource', ...
        '距离视图要求居中point源及direct；偏轴和激光暂不支持。');
    a=cfg.source.distance; b=cfg.observer.distance;
    scale=1/(1+b/a);
    step=2*deg2rad(cfg.observer.angleHalfRangeDeg)/(cfg.outputSize-1);
    kernels=coreObserverKernels(cfg,step); padding=kernels.paddingBins;
    reference=cfg;
    reference.observer.enabled=false;
    reference.source.type="plane";
    reference.source.wavefront="plane";
    if cfg.observer.referenceMethod=="fftBins"
        reference.receiverMethod="fftBins";
        reference.receiverFftPaddingFactor=cfg.observer.referencePaddingFactor;
    end
    reference.outputSize=cfg.outputSize+2*padding;
    reference.angleHalfRangeDeg=rad2deg(step*(reference.outputSize-1)/2)/scale;
    reference.metricAngleDeg=min(cfg.metricAngleDeg,reference.angleHalfRangeDeg);
    assert(reference.angleHalfRangeDeg<=10,'COEReflection:DistanceViewRange', ...
        '补边参考半角=%.4g°超过10°近轴范围；请缩小观察角域或源/核半径。', ...
        reference.angleHalfRangeDeg);
    % 不叠加球面波曲率；距离只通过几何映射进入，避免重复计入a。
    [mapped,baseline,heightEnabled]=coreSimulateHeightStates(reference);
    selected=padding+(1:cfg.outputSize);
    [energy,sourceEnergy,massError]=filterEnergy(mapped.spectralEnergy,kernels,selected);
    if heightEnabled
        [baselineEnergy,~,baselineMassError]=filterEnergy(baseline.spectralEnergy,kernels,selected);
    else
        baselineEnergy=energy; baselineMassError=massError;
    end
    view.model="Paraxial geometric scaling with finite-source angular averaging";
    view.assumptions="居中近轴；固定照明足迹；镜面贡献区域仅诊断；可选理想入瞳非相干角图PSF，无完整眼成像";
    view.sourceDistance=a;
    view.observerDistance=b;
    view.angularScale=scale;
    view.sourceDiameter=cfg.source.diameter; view.pupilDiameter=cfg.observer.pupilDiameter;
    view.screenScale=b*scale; % ab/(a+b)，屏幕投影长度/内禀角度。
    view.referenceCoordinate=mapped.coordinate(selected);
    view.coordinate=scale*view.referenceCoordinate;
    view.coordinateUnit="radian (observer angular coordinate, paraxial)";
    view.screenCoordinate=b*view.coordinate; % 视线与屏幕平面的交点，不是接收面光斑。
    view.pointSpectralEnergy=mapped.spectralEnergy(selected,selected,:);
    view.sourceAveragedSpectralEnergy=sourceEnergy;
    view.spectralEnergy=energy; view.baselineSpectralEnergy=baselineEnergy;
    weights=cfg.channelResponse.*cfg.spectralWeights(:);
    view.channelSignal=compose(energy,weights);
    view.baselineChannelSignal=compose(baselineEnergy,weights);
    view.binSolidAngle=step^2;
    view.spectralAngularDensity=energy/step^2;
    view.referenceThroughput=mapped.metrics.Throughput;
    view.capturedReferenceFraction=reshape(sum(sum(energy,1),2),[],1)./max(view.referenceThroughput,realmin);
    view.referenceConfig=reference;
    view.referenceMethod=cfg.observer.referenceMethod;
    view.config=cfg;
    view.heightPhaseEnabled=heightEnabled;
    view.referencePropagationRuns=mapped.propagationRuns;
    view.kernels=kernels; view.convolutionMassError=massError;
    view.baselineConvolutionMassError=baselineMassError;
    view.participation=coreObserverRegion(cfg,mapped);
    view.energyConvention="参考能量的角度推送；源概率核归一；Airy截尾核不补能量；非眼睛接收功率/W/m^2";
end

function signal=compose(energy,weights)
    shape=size(energy,1:2);
    signal=reshape(reshape(energy,[],size(weights,1))*weights,shape(1),shape(2),size(weights,2));
end

function [energy,sourceEnergy,errors]=filterEnergy(input,kernels,selected)
    energy=zeros(numel(selected),numel(selected),size(input,3)); sourceEnergy=energy;
    errors=zeros(size(input,3),2);
    for k=1:size(input,3)
        [source,errors(k,1)]=linearConvolve(input(:,:,k),kernels.sourceKernel);
        [pupil,errors(k,2)]=linearConvolve(source,kernels.pupilKernels{k});
        sourceEnergy(:,:,k)=source(selected,selected);
        energy(:,:,k)=pupil(selected,selected);
    end
end

function [same,error]=linearConvolve(input,kernel)
    if isscalar(kernel), same=input*kernel; error=0; return; end
    shape=size(input)+size(kernel)-1; padded=2.^nextpow2(shape);
    full=real(ifft2(fft2(input,padded(1),padded(2)).*fft2(kernel,padded(1),padded(2))));
    full=full(1:shape(1),1:shape(2));
    expected=sum(input,'all')*sum(kernel,'all'); error=abs(sum(full,'all')-expected);
    assert(error<=1e-11*max(1,expected),'COEReflection:ObserverConvolutionEnergy','角图线性卷积能量检查失败。');
    first=(size(kernel)-1)/2;
    same=max(0,full(first(1)+(1:size(input,1)),first(2)+(1:size(input,2))));
end
