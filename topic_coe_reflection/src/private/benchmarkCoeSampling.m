function report = benchmarkCoeSampling(cfg, stack, inputDir, outputDir, projectConfig)
%benchmarkCoeSampling 独立扫描采样参数；复用实际图片及统一相干膜堆。
% 光谱部分只计算区域系数/周期阶次，不构造数千波长的全物面传播数组。
% 空间部分使用生产传播入口；基准的平铺/抽样仅作用于内存副本。
    assert(stack.geometry.source=="images" && all(cfg.heightMap(:)==0) ...
        && isempty(cfg.complexCoefficient) && isinf(cfg.pupilRadius), ...
        'COEReflection:BenchmarkInput', '基准要求图片输入、零额外高度、空外部复系数及无限光瞳。');
    if ~isfolder(outputDir), mkdir(outputDir); end
    originalCfg=cfg;
    cfg.source=coeSourceDefaults();
    cfg.propagationDistance=100;
    cfg.wavelengths=(380:0.0625:780)*1e-9;
    fprintf('采样基准：区域光谱 0.0625 nm，D65 380–780 nm；空间扫描使用单色。\n');
    timer=tic;
    [model,geometry]=prepareCoeLayerStack(cfg,stack,inputDir);
    report.schemaVersion=1;
    report.createdAt=string(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
    report.matlabVersion=string(version);
    report.originalCfg=originalCfg; report.stack=stack; report.projectConfig=projectConfig;
    report.inputDir=string(inputDir);
    srcDir=fileparts(fileparts(mfilename('fullpath')));
    sourceFiles=[dir(fullfile(srcDir,'*.m'));dir(fullfile(srcDir,'private','*.m'))];
    report.sourceSnapshot=struct();
    for file=sourceFiles.'
        [~,name]=fileparts(file.name);
        report.sourceSnapshot.(name)=string(fileread(fullfile(file.folder,file.name)));
    end
    report.referenceWavelengthNm=cfg.wavelengths*1e9;
    report.regionTable=model.regionTable;
    report.coefficientLookup=model.coefficientLookup;
    report.materialTable=model.materialTable;
    report.protocol=struct('relativeTolerance',0.01,'absoluteTolerance',1e-4, ...
        'spectralCurveAbsoluteTolerance',0.01,'angleHalfRangeDeg',1, ...
        'coneHalfAngleDeg',0.8,'orderRoiHalfAnglesDeg',[0.08,0.10], ...
        'commonMapBins',200,'sourceDistanceM',0.1,'sourceDiameterM',0.2e-3, ...
        'note',"参照最细样本不等于精度保证；原生图片不支持更细真实几何的收敛声明。");
    cellSize=2*cfg.samplesPerPeriod;
    map=model.regionMap;
    periodic=all(mod(size(map),cellSize)==0) && ...
        isequal(map,repmat(map(1:cellSize,1:cellSize),size(map)/cellSize));
    report.twoPitchPeriodic=periodic;
    report.nativeRegionMap=map;
    report.spectral=spectralSweep(cfg,model,periodic,outputDir);
    writetable(report.spectral,fullfile(outputDir,'spectral.csv'));
    [report.ijpAnalyticError,report.ijpSpectral]=checkIjp(cfg,stack,inputDir,geometry);
    writetable(report.ijpSpectral,fullfile(outputDir,'ijp_spectral.csv'));
    report.spatial=table();
    report.commonMaps={};
    report.runConfigs={};
    warningState=warning('query','DisplayDiffraction:TruncatedEnergy');
    restoreWarning=onCleanup(@() warning(warningState));
    warning('off','DisplayDiffraction:TruncatedEnergy'); % 每行仍保留捕获率，不重复刷屏。
    base=cfg;
    base.angleHalfRangeDeg=1; base.metricAngleDeg=0.8;
    base.outputSize=401; base.receiverBinSamples=3;
    base.coeStack=model;
    base.coeStack.coefficientLookup=[]; base.coeStack.materialTable=table();
    % output：固定物面和 ±1°；bin：固定网格，扫描最少像素内积分点。
    for wavelength=[380,550,780]
        current=selectWavelength(base,model,wavelength);
        for count=[201,401,801,1601]
            current.outputSize=count;
            report=spatialRun(report,current,"output",count,outputDir);
        end
        % 同时检验日常 ±6°/801 设置；共同窗口仍固定为 ±1°便于对照。
        if deg2rad(6)*(1+1/400)<=current.wavelengths/(2*cfg.pitch/cfg.samplesPerPeriod)
            current.angleHalfRangeDeg=6;
            for count=[401,801,1201,1601]
                current.outputSize=count;
                report=spatialRun(report,current,"outputWide",count,outputDir);
            end
            current.angleHalfRangeDeg=1;
        end
        current.outputSize=401;
        for count=[1,3,5,7]
            current.receiverBinSamples=count;
            report=spatialRun(report,current,"receiverBin",count,outputDir);
        end
    end
    % 只对逐像素验证过的完整 2*Pitch 周期图样做重复数比较。
    if periodic
        unit=map(1:cellSize,1:cellSize);
        for wavelength=[380,550,780]
            for count=[8,16,32,64]
                current=selectWavelength(base,model,wavelength);
                current.outputSize=801;
                current.periods=[count,count];
                current.coeStack.regionMap=repmat(unit,count/2,count/2);
                report=spatialRun(report,current,"periods",count,outputDir);
            end
        end
    else
        fprintf('输入不是严格 2*Pitch 周期：跳过周期数扫描，避免制造缺陷超晶胞。\n');
    end
    % 原生图片的整数抽取（共同原点），不把插值图当成更精细版图。
    native=originalCfg.samplesPerPeriod;
    for wavelength=[380,550,780]
        for divisor=[4,2,1]
            if mod(native,divisor)~=0 || native/divisor<4, continue; end
            current=selectWavelength(base,model,wavelength);
            current.samplesPerPeriod=native/divisor;
            current.coeStack.regionMap=map(1:divisor:end,1:divisor:end);
            report=spatialRun(report,current,"objectRaster",native/divisor,outputDir);
        end
        % 零阶保持细化只检验同一硬像素图的数值积分，不能恢复亚像素轮廓。
        for factor=[1,2]
            current=selectWavelength(base,model,wavelength);
            current.samplesPerPeriod=native*factor;
            current.coeStack.regionMap=repelem(map,factor,factor);
            report=spatialRun(report,current,"pixelHold",native*factor,outputDir);
        end
    end
    % 有限源：固定 8*Pitch（非周期图使用原始全物面）、550 nm、201 网格。
    current=selectWavelength(base,model,550);
    current.outputSize=201;
    if periodic
        current.periods=[8,8];
        current.coeStack.regionMap=repmat(unit,4,4);
    end
    current.source.type="point";
    current.source.distance=report.protocol.sourceDistanceM;
    current.source.diameter=report.protocol.sourceDiameterM;
    for count=[3,5,9,17,25]
        current.source.samplesPerAxis=count;
        report=spatialRun(report,current,"source",count,outputDir);
    end
    report.spatial=compareSpatial(report.spatial,report.commonMaps,report.protocol);
    report.elapsedSeconds=toc(timer);
    writetable(report.spatial,fullfile(outputDir,'spatial.csv'));
    save(fullfile(outputDir,'sampling_convergence.mat'),'report','-v7.3');
    plotBenchmark(report,outputDir);
    fprintf('基准完成：%.1f s；%s\n',report.elapsedSeconds,outputDir);
end

function result = spectralSweep(cfg,model,periodic,outputDir)
    nm=cfg.wavelengths*1e9;
    orders=[0,0;1,0;0,1;1,1;2,0];
    regionWeights=model.regionTable.AreaFraction.';
    x=((0:size(model.regionMap,2)-1)-floor(size(model.regionMap,2)/2))/cfg.samplesPerPeriod;
    y=((0:size(model.regionMap,1)-1)-floor(size(model.regionMap,1)/2))/cfg.samplesPerPeriod;
    coefficients=zeros(size(orders,1),height(model.regionTable));
    for o=1:size(orders,1)
        phase=exp(-1i*pi*(orders(o,2)*y.'+orders(o,1)*x));
        for region=1:height(model.regionTable)
            coefficients(o,region)=sum(phase(model.regionMap==region-1))/numel(phase);
        end
    end
    reflectance=abs(model.coefficientLookup).^2;
    orderPower=abs(coefficients*model.coefficientLookup).^2;
    signals=[regionWeights*reflectance;orderPower];
    reference=integrateSpectrum(nm,signals);
    result=table();
    % 相移半格防止只在幸运的起点得出收敛结论；端点总是包含。
    for step=[5,2,1,0.5,0.25,0.125,0.0625]
        for shift=[0,0.5]
            if step==0.0625 && shift>0, continue; end
            selected=unique([380,380+shift*step:step:780,780]);
            index=round((selected-380)/0.0625)+1;
            coarse=signals(:,index);
            estimate=integrateSpectrum(selected,coarse);
            rCurve=interp1(selected,reflectance(:,index).',nm,'linear').';
            orderCurve=interp1(selected,orderPower(:,index).',nm,'linear').';
            xy=estimate.XYZ(:,1:2)./sum(estimate.XYZ,2);
            xyRef=reference.XYZ(:,1:2)./sum(reference.XYZ,2);
            active=reference.photopic>1e-8;
            row=table(step,shift,numel(selected), ...
                max(abs(rCurve-reflectance),[],'all'), ...
                max(abs(orderCurve-orderPower),[],'all'), ...
                estimate.radiant(1),estimate.photopic(1), ...
                abs(estimate.radiant(1)-reference.radiant(1)), ...
                abs(estimate.photopic(1)-reference.photopic(1)), ...
                max(abs(estimate.radiant(2:end)-reference.radiant(2:end))), ...
                max(abs(estimate.photopic(2:end)-reference.photopic(2:end))), ...
                max(vecnorm(xy(active,:)-xyRef(active,:),2,2)),periodic, ...
                'VariableNames',{'StepNm','GridShiftFraction','WavelengthCount', ...
                'MaxRegionRCurveAbsError','MaxOrderCurveAbsError','RadiantTotal','PhotopicTotal', ...
                'RadiantTotalAbsError','PhotopicTotalAbsError','MaxOrderRadiantAbsError', ...
                'MaxOrderPhotopicAbsError','MaxOrderXYDistance','TwoPitchPeriodic'});
            result=[result;row]; %#ok<AGROW>
        end
    end
    % 保存细网格原始曲线与阶次，作为可复核数据，不依赖显示增强。
    data=struct('wavelengthNm',nm,'regionReflectance',reflectance,'orders',orders, ...
        'orderPower',orderPower,'signals',signals,'integrated',reference, ...
        'interpretation',"阶次系数对应 2*Pitch 周期无限阵列；非周期输入仅为指定空间频率的诊断。");
    save(fullfile(outputDir,'spectral_reference.mat'),'data');
end

function values = integrateSpectrum(nm,signals)
    spectrum=loadCoeSpectrum(nm);
    values.radiant=signals*spectrum.energyWeights;
    values.photopic=signals*spectrum.photopicWeights;
    values.XYZ=signals*spectrum.xyzWeights;
end

function [error,rows] = checkIjp(cfg,stack,inputDir,geometry)
    error=NaN;
    rows=table();
    if ~any(stack.layerNames=="IJP"), return; end
    layer=stack.layers.IJP;
    if layer.method~="nk" || strlength(layer.nkFile)>0 || layer.n~=1.5 ...
            || layer.k~=0 || abs(layer.thickness-10e-6)>1e-15 ...
            || ~all(geometry.layerMasks.IJP,'all')
        return
    end
    stack.layerNames="IJP"; stack.layerLevels=1; stack.backgroundMode="medium";
    stack.incidentN=1; stack.exitN=1;
    film=prepareCoeLayerStack(cfg,stack,inputDir);
    phase=exp(4i*pi*1.5*10e-6./cfg.wavelengths);
    analytic=abs((-0.2+0.2*phase)./(1-0.04*phase)).^2;
    error=max(abs(abs(film.coefficientLookup).^2-analytic),[],'all');
    assert(error<1e-10,'COEReflection:BenchmarkAnalytic','IJP 单膜解析核对失败。');
    fprintf('10 μm / n=1.5 空气-膜-空气解析核对：%.3g。\n',error);
    nm=cfg.wavelengths*1e9;
    curve=abs(film.coefficientLookup(1,:)).^2;
    reference=integrateSpectrum(nm,curve);
    for step=[5,2,1,0.5,0.25,0.125]
        for shift=[0,0.5]
            selected=unique([380,380+shift*step:step:780,780]);
            index=round((selected-380)/0.0625)+1;
            sampled=curve(index);
            estimate=integrateSpectrum(selected,sampled);
            reconstructed=interp1(selected,sampled,nm,'linear');
            row=table(step,shift,max(abs(reconstructed-curve)), ...
                estimate.radiant,estimate.photopic, ...
                abs(estimate.radiant-reference.radiant), ...
                abs(estimate.photopic-reference.photopic), ...
                'VariableNames',{'StepNm','GridShiftFraction','MaxRCurveAbsError', ...
                'Radiant','Photopic','RadiantAbsError','PhotopicAbsError'});
            rows=[rows;row]; %#ok<AGROW>
        end
    end
end

function current = selectWavelength(base,model,nm)
    current=base;
    current.wavelengths=nm*1e-9;
    current.spectralWeights=1; current.channelResponse=1;
    [~,index]=min(abs(base.wavelengths*1e9-nm));
    current.coeStack.coefficientLookup=model.coefficientLookup(:,index);
end

function report = spatialRun(report,cfg,family,value,outputDir)
    fprintf('%s = %g，λ = %g nm\n',family,value,cfg.wavelengths*1e9);
    timer=tic;
    cfg.mask=ones(size(cfg.coeStack.regionMap)); % 图样已由统一膜堆读取，不调用历史几何分支。
    result=runCOEReflection(cfg);
    nm=cfg.wavelengths*1e9;
    energy=result.spectralEnergy;
    angle=rad2deg(result.coordinate);
    [xx,yy]=meshgrid(angle);
    center=rad2deg(cfg.wavelengths/(2*cfg.pitch));
    cone=sum(energy(hypot(xx,yy)<=0.8));
    order80=sum(energy(hypot(xx-center,yy)<=0.08));
    order100=sum(energy(hypot(xx-center,yy)<=0.10));
    axial=sum(energy(hypot(xx-2*center,yy)<=0.08));
    diagonal=sum(energy(hypot(xx-center,yy-center)<=0.08));
    step=angle(2)-angle(1);
    density=energy/deg2rad(step)^2;
    middle=(numel(angle)+1)/2;
    profile=density(middle,:);
    above=find(profile>=profile(middle)/2 & abs(angle)<0.2);
    width=NaN;
    if ~isempty(above), width=angle(above(end))-angle(above(1))+step; end
    targetEdges=linspace(-1,1,report.protocol.commonMapBins+1);
    left=angle-step/2; right=angle+step/2;
    overlap=max(0,min(targetEdges(2:end).',right)-max(targetEdges(1:end-1).',left))/step;
    common=overlap*energy*overlap.';
    row=table(family,value,nm,cfg.angleHalfRangeDeg,cfg.outputSize,cfg.receiverBinSamples, ...
        result.receiverSamplesPerBin,result.illumination.modeCount, ...
        cfg.periods(1),cfg.samplesPerPeriod,step, ...
        result.metrics.Throughput,sum(common,'all'),result.metrics.CapturedFraction, ...
        cone,order80,order100,axial,diagonal,max(density,[],'all'),width, ...
        result.metrics.ParsevalAbsoluteError,toc(timer), ...
        'VariableNames',{'Family','Value','WavelengthNm','HalfRangeDeg','OutputSize','MinimumBinSamples', ...
        'ActualBinSamples','SourceModeCount','PeriodsPerAxis','SamplesPerPitch','AngleStepDeg', ...
        'Total','CommonWindow','CapturedFraction','Cone80','OrderRoi80','OrderRoi100', ...
        'AxialOrderRoi80','DiagonalOrderRoi80','PeakPerRad2','CentralFwhmGridDeg','ParsevalAbsError','Seconds'});
    report.spatial=[report.spatial;row];
    report.commonMaps{end+1}=common;
    cfg.coeStack=rmfield(cfg.coeStack,{'materialTable','regionMap'});
    cfg.mask=[]; % 完整原图只保存一次；派生方式由 Family/Value 记录。
    report.runConfigs{end+1}=cfg;
    % 长扫描中保存已完成的行；中断后仍可审阅进度。
    writetable(report.spatial,fullfile(outputDir,'spatial_partial.csv'));
end

function rows = compareSpatial(rows,maps,protocol)
    rows.ReferenceValue=zeros(height(rows),1);
    rows.MaxMetricAbsError=zeros(height(rows),1);
    rows.CommonMapL1AbsError=zeros(height(rows),1);
    rows.IntegralWithinTolerance=false(height(rows),1);
    rows.CommonMapWithinTolerance=false(height(rows),1);
    rows.IsReference=false(height(rows),1);
    metrics={'Total','CommonWindow','Cone80','OrderRoi80','OrderRoi100','AxialOrderRoi80','DiagonalOrderRoi80'};
    for k=1:height(rows)
        group=find(rows.Family==rows.Family(k) & rows.WavelengthNm==rows.WavelengthNm(k));
        [~,index]=max(rows.Value(group)); reference=group(index);
        rows.ReferenceValue(k)=rows.Value(reference);
        rows.IsReference(k)=k==reference;
        delta=abs(rows{k,metrics}-rows{reference,metrics});
        rows.MaxMetricAbsError(k)=max(delta);
        rows.CommonMapL1AbsError(k)=sum(abs(maps{k}-maps{reference}),'all');
        rows.CommonMapWithinTolerance(k)=rows.CommonMapL1AbsError(k)<=protocol.absoluteTolerance ...
            +protocol.relativeTolerance*sum(maps{reference},'all');
        rows.IntegralWithinTolerance(k)=all(delta<=protocol.absoluteTolerance ...
            +protocol.relativeTolerance*abs(rows{reference,metrics}));
    end
end

function plotBenchmark(report,outputDir)
    fig=figure('Visible','off','Color','white','Position',[100,100,1400,650]);
    cleanup=onCleanup(@() close(fig));
    tiledlayout(2,4);
    nexttile;
    selected=report.spectral.StepNm>0.0625;
    semilogy(report.spectral.StepNm(selected),report.spectral.MaxRegionRCurveAbsError(selected),'o');
    xlabel('Spectral step / nm'); ylabel('Max regional R error'); grid on;
    families=["output","outputWide","receiverBin","periods","objectRaster","pixelHold","source"];
    for family=families
        nexttile; hold on;
        rows=report.spatial(report.spatial.Family==family,:);
        for nm=unique(rows.WavelengthNm).'
            chosen=rows(rows.WavelengthNm==nm,:);
            error=chosen.MaxMetricAbsError;
            error(chosen.IsReference)=NaN; % 参照样本自比较为零，不能伪装成收敛精度。
            line=semilogy(chosen.Value,error,'-o','DisplayName',string(nm)+' nm integral');
            mapError=chosen.CommonMapL1AbsError;
            mapError(chosen.IsReference)=NaN;
            semilogy(chosen.Value,mapError,'--','Color',line.Color, ...
                'DisplayName',string(nm)+' nm map L1');
        end
        title(family+" (ref "+string(max(rows.Value))+")");
        xlabel('Sampling parameter'); ylabel('Absolute difference');
        set(gca,'YScale','log');
        limits=ylim;
        decades=floor(log10(limits(1))):ceil(log10(limits(2)));
        ylim(10.^[decades(1),decades(end)]);
        yticks(10.^decades);
        grid on; legend('Location','best','FontSize',7);
    end
    exportgraphics(fig,fullfile(outputDir,'convergence.png'),'Resolution',150);
end
