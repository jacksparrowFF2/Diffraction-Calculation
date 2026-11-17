function report=benchDistanceScaling(cfg,stack,inputDir,outputDir,rgbOptions,projectConfig,saveResults)
%benchDistanceScaling 固定观察距离扫描光源，再固定光源距离扫描观察者。
% 复用真实图样和生产几何视图；量取550nm附近的标称Pitch第一阶候选峰。
    distances=cfg.observer.scanDistances;
    validateattributes(distances,{'numeric'},{'row','nonempty','real','finite','positive','increasing'});
    assert(numel(cfg.wavelengths)==3 && isequal(cfg.channelResponse,eye(3)), ...
        'COEReflection:DistanceBenchmarkSpectrum','距离扫描预览要求三个RGB代理通道。');
    [cfg.coeStack,geometry]=corePrepareLayerStack(cfg,stack,inputDir);
    cfg.mask=double(geometry.unionMask);
    count=numel(distances);
    a=[distances,repmat(cfg.source.distance,1,count)];
    b=[repmat(cfg.observer.distance,1,count),distances];
    views=cell(1,2*count); signals=views;
    predicted=zeros(2*count,1); measured=predicted; pointPeak=predicted; captured=predicted;
    comparable=false(2*count,1); diameter=predicted; periods=predicted;
    area=predicted; modeledArea=predicted; rayWeight=predicted; sourceWidth=predicted;
    [~,green]=min(abs(cfg.wavelengths-550e-9));
    for j=1:2*count
        fprintf('距离缩放 %d/%d：光源 %.1f cm，观察 %.1f cm\n',j,2*count,100*a(j),100*b(j));
        previous=find(a(1:j-1)==a(j) & b(1:j-1)==b(j),1);
        if isempty(previous)
            current=cfg; current.source.distance=a(j); current.observer.distance=b(j);
            views{j}=coreDistanceView(current);
        else
            views{j}=views{previous};
        end
        view=views{j};
        % 每单位角面积的物理参考信号先合并，再一次共同曝光。
        signals{j}=view.channelSignal/view.binSolidAngle;
        predicted(j)=rad2deg(view.angularScale*cfg.wavelengths(green)/cfg.pitch);
        coordinate=rad2deg(view.coordinate);
        candidates=find(coordinate>=0.7*predicted(j) & coordinate<=1.3*predicted(j));
        assert(~isempty(candidates),'COEReflection:DistanceBenchmarkWindow', ...
            '观察窗口未覆盖标称第一阶；请增大观察半角或提高outputSize。');
        centre=(numel(coordinate)+1)/2;
        profile=view.spectralEnergy(centre,candidates,green);
        assert(max(profile)>0,'COEReflection:DistanceBenchmarkPeak','候选范围内无信号。');
        [~,peak]=max(profile);
        measured(j)=coordinate(candidates(peak));
        [~,peak]=max(view.pointSpectralEnergy(centre,candidates,green));
        pointPeak(j)=coordinate(candidates(peak));
        sourceWidth(j)=rad2deg(view.kernels.sourceAngularDiameterRad);
        comparable(j)=sourceWidth(j)<0.8*predicted(j);
        if ~comparable(j), measured(j)=NaN; end % 源宽超过级次间距时不把混合峰伪报为单级峰。
        captured(j)=view.capturedReferenceFraction(green);
        diameter(j)=view.participation.diameter*1e3;
        periods(j)=view.participation.nominalPeriods;
        area(j)=view.participation.geometricArea*1e6;
        modeledArea(j)=view.participation.modeledArea*1e6;
        rayWeight(j)=view.participation.rayWeightInModel;
    end
    pooled=vertcat(signals{:});
    [rgb,~,mapping]=plotComposeReflectionRgb(pooled,pooled,rgbOptions,false);
    report=struct('model',views{1}.model,'assumptions',views{1}.assumptions, ...
        'baseConfig',cfg,'projectConfig',projectConfig,'layerGeometry',geometry, ...
        'views',{views},'mapping',mapping);
    report.peakNote="在0.7至1.3倍lambda/p预测位置内查找中心行峰；需输入确有该级次，不反推最小晶胞";
    report.baselineDisplayIndex=find(a==cfg.source.distance & b==cfg.observer.distance,1);
    report.summary=table(100*a(:),100*b(:),a(:)./(a(:)+b(:)),predicted,measured, ...
        measured-predicted,captured,'VariableNames', ...
        {'SourceDistanceCm','ObserverDistanceCm','AngularScale','PredictedFirstOrderDeg', ...
        'MeasuredPeakDeg','PeakErrorDeg','CapturedReferenceFraction'});
    report.summary.PointReferencePeakDeg=pointPeak;
    report.summary.PeakComparable=comparable;
    report.summary.SourceAngularDiameterDeg=sourceWidth;
    report.summary.PupilDiameterMm=repmat(cfg.observer.pupilDiameter*1e3,2*count,1);
    report.summary.GeometricDiameterMm=diameter;
    report.summary.NominalPeriods=periods;
    report.summary.GeometricAreaMm2=area;
    report.summary.ModeledAreaMm2=modeledArea;
    report.summary.RayWeightInModel=rayWeight;
    fig=figure('Name','COE 光源和观察距离缩放对比','Color','w','Position',[40,70,1500,760]);
    layout=tiledlayout(2,count,'TileSpacing','compact','Padding','compact');
    title(layout,'距离缩放观察角图 · 上排移动光源，下排移动观察者 · 全图共同曝光');
    subtitle(layout,sprintf('计算周期数 %d×%d Pitch；RGB代理；源直径%.1f mm，入瞳%.1f mm；%s', ...
        cfg.periods(1),cfg.periods(2),cfg.source.diameter*1e3,cfg.observer.pupilDiameter*1e3,cfg.observer.pupilModel));
    for j=1:2*count
        nexttile;
        rows=(j-1)*cfg.outputSize+(1:cfg.outputSize);
        coordinate=rad2deg(views{j}.coordinate);
        imagesc(coordinate,coordinate,rgb(rows,:,:)); axis image xy;
        xlabel('theta_x / °'); ylabel('theta_y / °');
        title(sprintf('a=%.0f cm, b=%.0f cm\ns=%.3f',100*a(j),100*b(j),views{j}.angularScale));
    end
    peakFigure=figure('Name','COE 距离缩放峰位验证','Color','w','Position',[80,100,1000,430]);
    peakLayout=tiledlayout(1,2,'TileSpacing','compact');
    title(peakLayout,sprintf('%d×%d Pitch · %.0f nm · +1级点源参考；有限源未分辨时不报告该级峰位', ...
        cfg.periods(1),cfg.periods(2),cfg.wavelengths(green)*1e9));
    for row=1:2
        nexttile; indices=(row-1)*count+(1:count);
        plot(100*distances,predicted(indices),'-','LineWidth',1.4); hold on;
        plot(100*distances,pointPeak(indices),'o','LineWidth',1.2);
        plot(100*distances,measured(indices),'x','LineWidth',1.2);
        if row==1, xlabel('光源距离 a / cm'); else, xlabel('观察距离 b / cm'); end
        ylabel('相对零级的观察角 / °'); grid on;
        legend('a/(a+b) · lambda/p','点源参考峰','最终图可比较峰','Location','best');
    end
    [regionFigure,stageFigure]=plotObserverDiagnostics(views,count,rgbOptions);
    if saveResults
        if ~isfolder(outputDir), mkdir(outputDir); end
        save(fullfile(outputDir,'distance_scaling_results.mat'),'report','-v7.3');
        writetable(report.summary,fullfile(outputDir,'distance_scaling_peaks.csv'));
        if isfield(projectConfig,'text')
            writelines(projectConfig.text,fullfile(outputDir,'config_snapshot.txt'),'Encoding','UTF-8');
        end
        exportgraphics(fig,fullfile(outputDir,'distance_scaling_comparison.png'),'Resolution',160);
        exportgraphics(peakFigure,fullfile(outputDir,'distance_scaling_peak_validation.png'),'Resolution',160);
        exportgraphics(regionFigure,fullfile(outputDir,'observer_participation_regions.png'),'Resolution',160);
        exportgraphics(stageFigure,fullfile(outputDir,'observer_source_pupil_stages.png'),'Resolution',160);
        writetable(report.summary,fullfile(outputDir,'observer_geometry_summary.csv'));
        if ~isempty(report.baselineDisplayIndex)
            rows=(report.baselineDisplayIndex-1)*cfg.outputSize+(1:cfg.outputSize);
            imwrite(flipud(rgb(rows,:,:)),fullfile(outputDir,'observer_baseline_rgb.png'));
        end
        fprintf('距离缩放结果：%s\n',outputDir);
    end
    disp(report.summary);
end
