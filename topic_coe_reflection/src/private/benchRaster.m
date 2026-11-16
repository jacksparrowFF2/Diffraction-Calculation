function report = benchRaster(cfg,stack,inputDir,outputDir,projectConfig)
%benchRaster 原生完整周期单元的物面采样收敛，不展开目标大阵列。
% 周期阶次用于初筛，有限物面生产传播用于复核；不改写输入图片。
    assert(stack.geometry.source=="images" && all(cfg.heightMap(:)==0) ...
        && isempty(cfg.complexCoefficient) && isinf(cfg.pupilRadius) && cfg.source.type=="plane", ...
        'COEReflection:BenchmarkInput','物面采样基准要求图片、平面波、零额外高度和无限光瞳。');
    if ~isfolder(outputDir), mkdir(outputDir); end
    started=tic;
    report.originalCfg=cfg; report.stack=stack; report.projectConfig=projectConfig;
    [cfg,native]=benchInputCell(cfg,stack,inputDir);
    periods=cfg.periods;
    report.protocol=struct('reference',"原生完整单元的无限周期阶次；不是连续几何真值", ...
        'relativeTolerance',0.01,'absoluteTolerance',1e-4,'mapL1Tolerance',0.01, ...
        'offsets',[0.25,0.5,0.75],'spectrumStepNm',5,'angleHalfRangeDeg',6, ...
        'metricAngleDeg',2.5,'requiredFinerPass',true);
    report.nativeSamplesPerPeriod=native; report.inputPeriods=periods;
    cfg.periods=periods; cfg.samplesPerPeriod=native;
    cfg.wavelengths=(380:5:780)*1e-9;
    cfg.spectralWeights=ones(size(cfg.wavelengths));
    cfg.channelResponse=ones(numel(cfg.wavelengths),1);
    stack.geometry.repeatToFit=false;
    stack.geometry.inputSamplesPerPeriod=native;
    [model,~]=corePrepareLayerStack(cfg,stack,inputDir);
    nativeMap=model.regionMap;
    spectrum=coreLoadSpectrum(cfg.wavelengths*1e9);
    % 所有候选在同一组物理阶次上比较，频率由真实单元宽度决定。
    maxOrder=ceil(deg2rad(6).*(periods*cfg.pitch)/min(cfg.wavelengths));
    [mx,my]=meshgrid(-maxOrder(1):maxOrder(1),-maxOrder(2):maxOrder(2));
    orders=[mx(:),my(:)];
    [reference,referencePower]=measureRaster(nativeMap,model,orders,cfg,periods);
    report.wavelengthNm=cfg.wavelengths*1e9; report.orders=orders;
    report.reference=reference; report.regionTable=model.regionTable;
    report.coefficientLookup=model.coefficientLookup;
    report.rows=table(); report.spectral=struct('samples',{},'offset',{},'values',{});
    counts=unique([40,50,64,80,100,125,160,200,native]);
    counts=counts(counts<=native & counts>2*max(maxOrder./periods));
    for count=counts
        for offset=report.protocol.offsets
            rows=floor(((0:periods(2)*count-1)+offset)*native/count)+1;
            cols=floor(((0:periods(1)*count-1)+offset)*native/count)+1;
            map=nativeMap(rows,cols);
            [values,power]=measureRaster(map,model,orders,cfg,periods);
            error=abs(values-reference);
            allowed=1e-4+0.01*abs(reference);
            l1=max(sum(abs(power-referencePower),1)./max(sum(referencePower,1),1e-12));
            pass=all(error<=allowed,'all') && l1<=0.01;
            target=report.originalCfg.periods*count;
            row=table(count,offset,max(error,[],'all'),max(error./allowed,[],'all'),l1, ...
                max(abs((values-reference)*spectrum.energyWeights)), ...
                max(abs((values-reference)*spectrum.photopicWeights)),pass, ...
                prod(target)*16/2^30,prod(target)*4/2^30, ...
                'VariableNames',{'SamplesPerPitch','PixelOffset','MaxCurveAbsError', ...
                'MaxToleranceRatio','AngularL1RelativeError','D65RadiantAbsError', ...
                'D65PhotopicAbsError','Pass','OneComplexGridGiB','RegionMapGiB'});
            report.rows=[report.rows;row];
            report.spectral(end+1)=struct('samples',count,'offset',offset,'values',values);
        end
        fprintf('物面采样：%g / %g pixels/Pitch，%.1f s\n',count,native,toc(started));
    end
    report.periodicCandidate=native;
    for count=counts(counts<native)
        finer=report.rows.SamplesPerPitch>=count;
        if all(report.rows.Pass(finer))
            report.periodicCandidate=count; break
        end
    end
    writetable(report.rows,fullfile(outputDir,'raster_summary.csv'));
    save(fullfile(outputDir,'raster_periodic_screen.mat'),'report','-v7.3');
    % 完整原图单元上的有限阵列复核，覆盖全谱；与周期初筛分开记录。
    report.production=validateProduction(cfg,stack,inputDir,report.periodicCandidate,native,outputDir);
    report.recommendedSamplesPerPeriod=native;
    if report.periodicCandidate<native && all(report.production.Pass)
        report.recommendedSamplesPerPeriod=report.periodicCandidate;
    end
    report.recommendationScope="平面波、当前材料/厚度、380:5:780nm、±6°；有限物面复核仅原图4Pitch，64Pitch窄峰形状仍需另验";
    report.elapsedSeconds=toc(started);
    writetable(report.rows,fullfile(outputDir,'raster_summary.csv'));
    writetable(report.production,fullfile(outputDir,'production_validation.csv'));
    save(fullfile(outputDir,'raster_convergence.mat'),'report','-v7.3');
    f=figure('Visible','off'); cleanup=onCleanup(@() close(f));
    tiledlayout(1,2);
    nexttile; hold on;
    for offset=report.protocol.offsets
        selected=report.rows.PixelOffset==offset;
        plot(report.rows.SamplesPerPitch(selected),report.rows.MaxToleranceRatio(selected),'-o', ...
            'DisplayName',sprintf('pixel offset %.2f',offset));
    end
    yline(1,'--'); xlabel('pixels / Pitch'); ylabel('error / tolerance'); legend('Location','best'); grid on;
    nexttile; selected=report.rows.PixelOffset==0.5;
    plot(report.rows.SamplesPerPitch(selected),report.rows.OneComplexGridGiB(selected),'-o');
    xlabel('pixels / Pitch'); ylabel('one complex grid / GiB'); grid on;
    exportgraphics(f,fullfile(outputDir,'raster_convergence.png'),'Resolution',150);
    fprintf('物面采样建议：%g pixels/Pitch；%s\n',report.recommendedSamplesPerPeriod,outputDir);
end

function [metrics,power]=measureRaster(map,model,orders,cfg,periods)
    shape=size(map);
    assert(all(max(abs(orders),[],1)<shape([2,1])/2), ...
        'COEReflection:ReceiverAliasing','候选网格不能覆盖±6°所有阶次。');
    coefficients=complex(zeros(size(orders,1),height(model.regionTable)));
    index=sub2ind(shape,mod(orders(:,2),shape(1))+1,mod(orders(:,1),shape(2))+1);
    fractions=zeros(1,height(model.regionTable));
    for region=1:height(model.regionTable)
        mask=map==region-1;
        transform=fft2(double(mask))/numel(map);
        coefficients(:,region)=transform(index);
        fractions(region)=mean(mask,'all');
    end
    power=abs(coefficients*model.coefficientLookup).^2;
    ax=orders(:,1)./(periods(1)*cfg.pitch)*cfg.wavelengths;
    ay=orders(:,2)./(periods(2)*cfg.pitch)*cfg.wavelengths;
    window=abs(ax)<=deg2rad(6) & abs(ay)<=deg2rad(6);
    radius=hypot(ax,ay);
    power(~window)=0;
    metrics=[fractions*abs(model.coefficientLookup).^2;sum(power,1); ...
        sum(power.*(radius<=deg2rad(2.5)),1);sum(power.*(radius>deg2rad(0.3)),1)];
end

function rows=validateProduction(cfg,stack,inputDir,candidate,native,outputDir)
    % 逐波长流式验证，避免存储完整物面×波长；输出共同±6°角度能量网格。
    rows=table();
    cfg.outputSize=121; cfg.angleHalfRangeDeg=6; cfg.metricAngleDeg=2.5;
    cfg.channelResponse=1; cfg.spectralWeights=1;
    warningState=warning('query','DisplayDiffraction:TruncatedEnergy');
    cleanup=onCleanup(@() warning(warningState));
    warning('off','DisplayDiffraction:TruncatedEnergy');
    wavelengths=380:5:780;
    if candidate==native, wavelengths=550; end
    for nm=wavelengths
        cfg.wavelengths=nm*1e-9;
        cfg.samplesPerPeriod=native;
        reference=runCOEReflection(cfg,stack,inputDir);
        cfg.samplesPerPeriod=candidate;
        value=runCOEReflection(cfg,stack,inputDir);
        [xx,yy]=meshgrid(value.coordinate);
        roi=hypot(xx,yy)<=deg2rad(2.5);
        a=[reference.metrics.Throughput;sum(reference.spectralEnergy,'all');sum(reference.spectralEnergy(roi))];
        b=[value.metrics.Throughput;sum(value.spectralEnergy,'all');sum(value.spectralEnergy(roi))];
        ratio=max(abs(b-a)./(1e-4+0.01*abs(a)));
        l1=sum(abs(value.spectralEnergy-reference.spectralEnergy),'all')/max(sum(reference.spectralEnergy,'all'),1e-12);
        rows=[rows;table(nm,ratio,l1,ratio<=1 && l1<=0.01, ...
            max(reference.metrics.ParsevalAbsoluteError),max(value.metrics.ParsevalAbsoluteError), ...
            'VariableNames',{'WavelengthNm','MaxToleranceRatio','AngularL1RelativeError','Pass', ...
            'ReferenceParsevalError','CandidateParsevalError'})]; %#ok<AGROW>
        if mod(nm-380,50)==0
            writetable(rows,fullfile(outputDir,'production_validation_partial.csv'));
            fprintf('有限物面生产复核：%g nm\n',nm);
        end
    end
end
