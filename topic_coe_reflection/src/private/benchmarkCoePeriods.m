function report = benchmarkCoePeriods(cfg,stack,inputDir,outputDir,projectConfig)
%benchmarkCoePeriods 已验证完整晶胞的有限重复阵列效率，逐波长与无限周期比较。
% 用晶胞非周期自相关与有限阵列三角权重精确积分角度矩形；不展开大物面。
% 与生产代码采用相同的离散物面点求积，不增加像素孔径或删除膜堆相位。
    assert(stack.geometry.source=="images" && all(cfg.heightMap(:)==0) ...
        && isempty(cfg.complexCoefficient) && isinf(cfg.pupilRadius), ...
        'COEReflection:BenchmarkInput','Pitch基准要求图片输入、零额外高度、无外部复系数及无限光瞳。');
    if ~isfolder(outputDir), mkdir(outputDir); end
    started=tic;
    originalCfg=cfg;
    cfg.wavelengths=(380:5:780)*1e-9;
    cfg.spectralWeights=ones(size(cfg.wavelengths));
    cfg.channelResponse=ones(numel(cfg.wavelengths),1);
    cfg.source=coeSourceDefaults(); % 平面波隔离重复数效应；不将有限源改成平面波结论。
    cfg.propagationDistance=1e4; % 只作角分布参考；效率与该距离无关。
    cfg.angleHalfRangeDeg=6; cfg.outputSize=1201;
    [model,~]=prepareCoeLayerStack(cfg,stack,inputDir);
    q=2*cfg.samplesPerPeriod;
    map=model.regionMap;
    assert(all(mod(size(map),q)==0) && ...
        isequal(map,repmat(map(1:q,1:q),size(map)/q)), ...
        'COEReflection:BenchmarkPeriodicity','图片不是严格2×Pitch周期，不能据此推断无限周期收敛。');
    assert(deg2rad(6)*(1+1/1200)<=min(cfg.wavelengths)/(2*cfg.pitch/cfg.samplesPerPeriod), ...
        'COEReflection:ReceiverAliasing','±6°超过物面采样带宽。');
    unit=map(1:q,1:q);
    counts=[8,16,24,32,48,64,96,128,192,256];
    orders=[0,0;1,0;0,1;1,1;2,0;0,2;2,2;4,0;4,4;8,0;12,0];
    halves=[0.08,0.10];
    metricNames="window6";
    for half=halves
        for k=1:size(orders,1)
            metricNames(end+1)=sprintf('m%d_n%d_h%03d',orders(k,1),orders(k,2),round(half*1000)); %#ok<AGROW>
        end
    end
    report.schemaVersion=1;
    report.createdAt=string(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
    report.matlabVersion=string(version);
    report.originalCfg=originalCfg; report.stack=stack; report.projectConfig=projectConfig;
    report.inputDir=string(inputDir); report.nativeRegionMap=map; report.unitRegionMap=unit;
    report.regionTable=model.regionTable; report.materialTable=model.materialTable;
    report.coefficientLookup=model.coefficientLookup;
    report.countsPitch=counts; report.orders=orders; report.roiHalfAnglesDeg=halves;
    report.wavelengthNm=cfg.wavelengths*1e9; report.metricNames=metricNames;
    report.protocol=struct('sourceType',"plane",'spectrum',"D65",'wavelengthStepNm',5, ...
        'outputHalfAngleDeg',6,'outputSize',1201,'relativeTolerance',0.01, ...
        'absoluteTolerance',1e-4,'requiredConsecutive',2, ...
        'roiShape',"square",'reference',"infinite periodic diffraction-order efficiency", ...
        'note',"每轴Pitch计数；不以峰高或点对点峰形作为周期效率收敛判据。");
    report.efficiency=zeros(numel(counts),numel(cfg.wavelengths),numel(metricNames));
    report.infiniteEfficiency=zeros(numel(cfg.wavelengths),numel(metricNames));
    report.throughput=zeros(numel(cfg.wavelengths),1);
    report.parsevalError=zeros(numel(counts),numel(cfg.wavelengths));
    dx=cfg.pitch/cfg.samplesPerPeriod;
    lag=(-(q-1):(q-1)).';
    % 分别建立x、y方向的晶胞内相位列。
    phaseX=exp(-2i*pi*lag*orders(:,1).'/q);
    phaseY=exp(-2i*pi*lag*orders(:,2).'/q);
    sourceFiles=[dir(fullfile(fileparts(fileparts(mfilename('fullpath'))),'*.m')); ...
        dir(fullfile(fileparts(mfilename('fullpath')),'*.m'))];
    for file=sourceFiles.'
        [~,name]=fileparts(file.name);
        report.sourceSnapshot.(name)=string(fileread(fullfile(file.folder,file.name)));
    end
    for w=1:numel(cfg.wavelengths)
        lambda=cfg.wavelengths(w);
        values=model.coefficientLookup(:,w);
        field=reshape(values(double(unit)+1),size(unit));
        correlation=unitCorrelation(field);
        power=abs(fftshift(fft2(field))/q^2).^2;
        frequencies=(-q/2:q/2-1)/(q*dx);
        inWindow=abs(lambda*frequencies)<=deg2rad(6);
        report.infiniteEfficiency(w,1)=sum(power(inWindow,inWindow),'all');
        report.throughput(w)=mean(abs(field).^2,'all');
        for h=1:numel(halves)
            for o=1:size(orders,1)
                center=rad2deg(lambda*orders(o,:)/(q*dx));
                assert(all(abs(center)+halves(h)<6), ...
                    'COEReflection:BenchmarkRoi','阶次ROI越出±6°。');
                index=1+(h-1)*size(orders,1)+o;
                report.infiniteEfficiency(w,index)=power(q/2+1+orders(o,2),q/2+1+orders(o,1));
            end
        end
        for n=1:numel(counts)
            repetitions=counts(n)/2;
            area=(counts(n)*cfg.pitch)^2;
            scale=dx^4/area;
            kernel=integralKernel(q,repetitions,dx,deg2rad(6)/lambda,0);
            report.efficiency(n,w,1)=real(kernel.'*correlation*kernel)*scale;
            for h=1:numel(halves)
                base=integralKernel(q,repetitions,dx,deg2rad(halves(h))/lambda,0);
                kx=base.*phaseX; ky=base.*phaseY;
                integrated=real(sum(ky.*(correlation*kx),1))*scale;
                selected=1+(h-1)*size(orders,1)+(1:size(orders,1));
                report.efficiency(n,w,selected)=integrated;
            end
            nyquist=integralKernel(q,repetitions,dx,1/(2*dx),0);
            total=real(nyquist.'*correlation*nyquist)*scale;
            report.parsevalError(n,w)=abs(total-report.throughput(w));
        end
        if mod(w,10)==1, fprintf('Pitch积分：%g nm，8–256 Pitch\n',lambda*1e9); end
    end
    assert(max(report.parsevalError,[],'all')<1e-9,'COEReflection:BenchmarkParseval','周期积分Parseval检查失败。');
    assert(min(report.efficiency,[],'all')>=-1e-10,'COEReflection:BenchmarkPassivity','积分出现负能量。');
    spectrum=loadCoeSpectrum(report.wavelengthNm);
    report.radiant=reshape(permute(report.efficiency,[1,3,2]),[],numel(cfg.wavelengths))*spectrum.energyWeights;
    report.radiant=reshape(report.radiant,numel(counts),[]);
    report.photopic=reshape(permute(report.efficiency,[1,3,2]),[],numel(cfg.wavelengths))*spectrum.photopicWeights;
    report.photopic=reshape(report.photopic,numel(counts),[]);
    report.infiniteRadiant=spectrum.energyWeights.'*report.infiniteEfficiency;
    report.infinitePhotopic=spectrum.photopicWeights.'*report.infiniteEfficiency;
    report.summary=summarizePeriods(report);
    writetable(report.summary,fullfile(outputDir,'pitch_summary.csv'));
    report.validation=validatePeriodicIntegral(cfg,model,unit,dx);
    writetable(report.validation,fullfile(outputDir,'production_validation.csv'));
    report.preview=periodPreviews(cfg,model,unit,dx,report.summary,outputDir);
    report.elapsedSeconds=toc(started);
    save(fullfile(outputDir,'pitch_convergence.mat'),'report','-v7.3');
    plotPeriods(report,outputDir);
    fprintf('Pitch基准完成：%.1f s；%s\n',report.elapsedSeconds,outputDir);
    disp(report.summary);
end

function correlation = unitCorrelation(field)
    q=size(field,1);
    padded=fft2(field,2*q,2*q);
    correlation=fftshift(ifft2(abs(padded).^2));
    correlation=correlation(2:end,2:end); % lags=-(q-1):(q-1)，不重复端点。
end

function kernel = integralKernel(q,repetitions,dx,halfWidth,center)
% 对s个晶胞内部像素差，累加各晶胞对的三角权重。
% integral exp(-2*pi*i*f*position) df = 2*h*sinc(2*h*position)*exp(-2*pi*i*c*position)。
    lag=(-(q-1):(q-1)).';
    shift=-(repetitions-1):(repetitions-1);
    position=dx*(lag+q*shift);
    argument=2*halfWidth*position;
    sincValue=ones(size(argument));
    active=abs(argument)>1e-14;
    sincValue(active)=sin(pi*argument(active))./(pi*argument(active));
    weights=(repetitions-abs(shift)).*2*halfWidth.*sincValue;
    kernel=weights*exp(-2i*pi*position(1,:).'*center); % 下面补足各lag相对于首行的相位。
    kernel=kernel.*exp(-2i*pi*(lag-lag(1))*dx*center);
end

function energy = periodicBins(field,countPitch,dx,lambda,angleDeg)
% 完整输出bin的矩形精确积分，避免大阵列的显式细网格。
    q=size(field,1);
    coordinate=deg2rad(angleDeg);
    step=coordinate(2)-coordinate(1);
    kernel=integralKernel(q,countPitch/2,dx,step/(2*lambda),coordinate/lambda);
    correlation=unitCorrelation(field);
    energy=real(kernel.'*correlation*kernel)*dx^4/(countPitch*q/2*dx)^2;
end

function summary = summarizePeriods(report)
    summary=table();
    reference=report.infiniteEfficiency;
    for n=1:numel(report.countsPitch)
        values=reshape(report.efficiency(n,:,:),size(reference));
        delta=abs(values-reference);
        tolerance=report.protocol.absoluteTolerance+report.protocol.relativeTolerance*abs(reference);
        radiantDelta=abs(report.radiant(n,:)-report.infiniteRadiant);
        photopicDelta=abs(report.photopic(n,:)-report.infinitePhotopic);
        broadbandPassed=all(radiantDelta<=report.protocol.absoluteTolerance ...
            +report.protocol.relativeTolerance*abs(report.infiniteRadiant)) ...
            && all(photopicDelta<=report.protocol.absoluteTolerance ...
            +report.protocol.relativeTolerance*abs(report.infinitePhotopic));
        [worst,index]=max(delta./tolerance,[],'all','linear');
        [w,m]=ind2sub(size(delta),index);
        row=table(report.countsPitch(n),report.countsPitch(n)*report.originalCfg.pitch*1e3, ...
            max(delta,[],'all'),max(radiantDelta),max(photopicDelta),worst, ...
            report.wavelengthNm(w),report.metricNames(m),all(delta<=tolerance,'all'),broadbandPassed, ...
            'VariableNames',{'PitchPerAxis','WidthMm','MaxSpectralAbsError', ...
            'MaxRadiantAbsError','MaxPhotopicAbsError','WorstToleranceRatio', ...
            'WorstWavelengthNm','WorstMetric','AllWavelengthsPassed','BroadbandPassed'});
        summary=[summary;row]; %#ok<AGROW>
    end
    % 不能凭一个偶然过关的样本推荐：要求该点起所有更大样本通过且至少再有一个样本。
    summary.Recommended=false(height(summary),1);
    summary.BroadbandRecommended=false(height(summary),1);
    for n=1:height(summary)-1
        if all(summary.AllWavelengthsPassed(n:end)) && all(summary.BroadbandPassed(n:end))
            summary.Recommended(n)=true;
            break
        end
    end
    for n=1:height(summary)-1
        if all(summary.BroadbandPassed(n:end))
            summary.BroadbandRecommended(n)=true;
            break
        end
    end
end

function validation = validatePeriodicIntegral(cfg,model,unit,dx)
    validation=table();
    state=warning('query','DisplayDiffraction:TruncatedEnergy');
    cleanup=onCleanup(@() warning(state));
    warning('off','DisplayDiffraction:TruncatedEnergy');
    for nm=[380,550,780]
        [~,index]=min(abs(cfg.wavelengths*1e9-nm));
        values=model.coefficientLookup(:,index);
        field=reshape(values(double(unit)+1),size(unit));
        for count=[8,16]
            current=cfg; current.wavelengths=nm*1e-9;
            current.spectralWeights=1; current.channelResponse=1;
            current.periods=[count,count]; current.angleHalfRangeDeg=1;
            current.metricAngleDeg=0.8; current.outputSize=201; current.receiverBinSamples=9;
            current.coeStack=model; current.coeStack.regionMap=repmat(unit,count/2,count/2);
            current.coeStack.coefficientLookup=values;
            current.mask=ones(size(current.coeStack.regionMap));
            production=runCOEReflection(current);
            exact=periodicBins(field,count,dx,nm*1e-9,rad2deg(production.coordinate));
            difference=sum(abs(exact-production.spectralEnergy),'all');
            relative=difference/sum(exact,'all');
            assert(relative<1e-3,'COEReflection:BenchmarkValidation','周期积分与生产传播不符。');
            row=table(count,nm,relative,max(abs(exact-production.spectralEnergy),[],'all'), ...
                'VariableNames',{'PitchPerAxis','WavelengthNm','MapL1RelativeError','MaxBinAbsError'});
            validation=[validation;row]; %#ok<AGROW>
        end
    end
end

function preview = periodPreviews(cfg,model,unit,dx,summary,outputDir)
    candidate=summary.PitchPerAxis(summary.Recommended);
    if isempty(candidate), candidate=summary.PitchPerAxis(end); end
    broadband=summary.PitchPerAxis(summary.BroadbandRecommended);
    if isempty(broadband), broadband=candidate; end
    counts=unique([16,64,broadband,candidate]);
    [~,index]=min(abs(cfg.wavelengths*1e9-550));
    values=model.coefficientLookup(:,index);
    field=reshape(values(double(unit)+1),size(unit));
    angles=linspace(-6,6,1201);
    preview=struct('wavelengthNm',550,'pitchPerAxis',counts,'angleDeg',angles, ...
        'energy',{{}},'sourceType',"plane",'referenceDensity',0);
    for count=counts
        fprintf('Pitch角图：%d Pitch，550 nm，±6°/1201\n',count);
        preview.energy{end+1}=periodicBins(field,count,dx,550e-9,angles);
    end
    step=deg2rad(angles(2)-angles(1));
    preview.referenceDensity=max(cellfun(@(x) max(x,[],'all'),preview.energy),[],'all')/step^2;
    fig=figure('Visible','off','Color','white','Position',[100,100,1400,450]);
    cleanup=onCleanup(@() close(fig));
    layout=tiledlayout(fig,1,numel(counts),'Padding','loose','TileSpacing','compact');
    for k=1:numel(counts)
        ax=nexttile(layout);
        db=10*log10(max(preview.energy{k}/step^2/preview.referenceDensity,1e-5));
        imagesc(ax,angles,angles,db); axis(ax,'image'); set(ax,'YDir','normal');
        xlim(ax,[-1,1]); ylim(ax,[-1,1]); clim(ax,[-50,0]); colormap(ax,gray(256));
        title(ax,sprintf('%d x %d Pitch',counts(k),counts(k)));
        xlabel(ax,'theta_x / deg'); ylabel(ax,'theta_y / deg');
    end
    cb=colorbar(ax); cb.Layout.Tile='east'; cb.Label.String='Common energy-density reference / dB';
    title(layout,'550 nm | plane wave | +/-6 deg / 1201 | view +/-1 deg');
    drawnow;
    exportgraphics(layout,fullfile(outputDir,'pitch_images_550nm.png'),'Resolution',150);
end

function plotPeriods(report,outputDir)
    fig=figure('Visible','off','Color','white','Position',[100,100,1150,450]);
    cleanup=onCleanup(@() close(fig));
    layout=tiledlayout(fig,1,2,'Padding','loose');
    ax=nexttile(layout);
    semilogy(ax,report.summary.PitchPerAxis,report.summary.WorstToleranceRatio,'-o');
    yline(ax,1,'--','Tolerance'); grid(ax,'on');
    xlabel(ax,'Pitch per axis'); ylabel(ax,'Worst absolute difference / tolerance');
    title(ax,'All 81 wavelengths and fixed order ROIs');
    ax=nexttile(layout);
    semilogy(ax,report.summary.PitchPerAxis,report.summary.MaxRadiantAbsError,'-o', ...
        report.summary.PitchPerAxis,report.summary.MaxPhotopicAbsError,'-s');
    grid(ax,'on'); legend(ax,'D65 radiant','D65 photopic');
    xlabel(ax,'Pitch per axis'); ylabel(ax,'Max absolute reflectance difference');
    title(ax,'Relative to the infinite periodic array');
    drawnow;
    exportgraphics(layout,fullfile(outputDir,'pitch_convergence.png'),'Resolution',150);
end
