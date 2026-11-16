function report = benchPlaneWaveComparison(cfg, stack, inputDir, outputDir, rgbOptions, projectConfig)
%benchPlaneWaveComparison 真实完整图片单元与显式有限周期阵列的平面波RGB对比。
% 不裁切、缩放或改写原图；所有图采用共同角网格、共同曝光及单位入射功率归一化。
    assert(cfg.source.type=="plane" && all(cfg.heightMap(:)==0) ...
        && isempty(cfg.complexCoefficient) && isinf(cfg.pupilRadius), ...
        'COEReflection:ComparisonInput','要求平面波、零额外高度、无限光瞳及图层复系数。');
    assert(isequal(string(stack.layerNames),["Anode","BM"]), ...
        'COEReflection:ComparisonLayers','本结构对比仅使用Anode/BM图层。');
    assert(isequal(cfg.channelResponse,eye(3)) && numel(cfg.wavelengths)==3 ...
        && all(diff(cfg.wavelengths)<0), 'COEReflection:ComparisonRGB', ...
        '要求按R/G/B排列的三个代理波长和eye(3)通道响应。');
    [localCfg,native] = benchInputCell(cfg,stack,inputDir);
    assert(cfg.samplesPerPeriod==native && any(cfg.periods>localCfg.periods), ...
        'COEReflection:ComparisonSampling','本对比保留原生采样，阵列须大于原图单元。');
    assert(stack.geometry.repeatToFit, 'COEReflection:ComparisonRepeat', ...
        '阵列比较必须显式开启repeatToFit。');
    if ~isfolder(outputDir), mkdir(outputDir); end
    report.protocol = "平面波；原图完整单元与有限平铺阵列；三波长RGB代理；共同曝光；无相机模型";
    report.projectConfig=projectConfig; report.stack=stack;
    report.inputDir=string(inputDir);
    cases={localCfg,cfg}; labels=["local","array"];
    started=tic;
    summary=table();
    for j=1:2
        current=cases{j};
        fprintf('平面波对比 %s：%g×%g Pitch，%g点/Pitch\n', ...
            labels(j),current.periods,current.samplesPerPeriod);
        value=runCOEReflection(current,stack,inputDir);
        assert(value.propagationRuns==1 && all(value.metrics.ParsevalAbsoluteError<1e-10), ...
            'COEReflection:ComparisonEnergy','单状态传播或Parseval检查失败。');
        if j==1
            localMask=value.layerGeometry.layerMasks.Anode;
            % 保存实际有效复反射区域；BM在上层为零反射边界。
            model=value.config.coeStack;
            localEffective=reshape(model.coefficientLookup(double(model.regionMap)+1,2),size(model.regionMap));
            clear model
        end
        data=struct('cfg',current,'coordinate',value.coordinate, ...
            'spectralEnergy',value.spectralEnergy,'channelSignal',value.channelSignal, ...
            'metrics',value.metrics,'incidentArea',value.incidentArea, ...
            'incidentPower',value.incidentPower,'repeatCounts',value.layerGeometry.repeatCounts, ...
            'illuminationAperture',value.illuminationAperture, ...
            'inputPeriods',value.layerGeometry.inputPeriods, ...
            'receiverSamplesPerBin',value.receiverSamplesPerBin);
        report.(labels(j))=data;
        row=value.metrics;
        row.Case=repmat(labels(j),height(row),1);
        row.PeriodsX=repmat(current.periods(1),height(row),1);
        row.PeriodsY=repmat(current.periods(2),height(row),1);
        row.IncidentAreaM2=repmat(value.incidentArea,height(row),1);
        summary=[summary;row]; %#ok<AGROW>
        clear value
        fprintf('完成 %s，累计 %.1f s\n',labels(j),toc(started));
    end
    assert(isequal(report.local.coordinate,report.array.coordinate), ...
        'COEReflection:ComparisonGrid','两组角网格不一致。');
    assert(report.local.illuminationAperture.coverageFraction<1 ...
        || report.array.illuminationAperture.coverageFraction<1 ...
        || max(abs(report.local.metrics.Throughput-report.array.metrics.Throughput))<1e-12, ...
        'COEReflection:ComparisonTiling','完整平铺后的反射率不一致。');
    [localRGB,arrayRGB,report.mapping]=plotComposeReflectionRgb( ...
        report.local.channelSignal,report.array.channelSignal,rgbOptions,true);
    report.local.rgb=localRGB; report.array.rgb=arrayRGB;
    report.summary=summary; report.elapsedSeconds=toc(started);
    report.local.anodeMask=logical(localMask);
    report.local.effectiveReflection=localEffective;
    writetable(summary,fullfile(outputDir,'plane_wave_summary.csv'));
    imwrite(flipud(localRGB),fullfile(outputDir,'local_rgb.png'));
    imwrite(flipud(arrayRGB),fullfile(outputDir,'array_rgb.png'));
    save(fullfile(outputDir,'plane_wave_comparison.mat'),'report','-v7.3');
    renderComparison(report,outputDir);
    fprintf('平面波对比完成：%s\n',outputDir);
end

function renderComparison(report,outputDir)
    coord=rad2deg(report.local.coordinate);
    half=report.local.cfg.angleHalfRangeDeg;
    zoomHalf=min(1.5,half);
    cases={report.local,report.array};
    labels=["完整原图局部单元","周期阵列"];
    % 单波长强度参考共用，不分别归一化；RGB采用同一个映射参考值。
    greenReference=max([report.local.spectralEnergy(:,:,2),report.array.spectralEnergy(:,:,2)],[],'all');
    fig=figure('Color','w','Position',[60,60,1500,820], ...
        'Name','Anode/BM 平面波：局部单元与周期阵列');
    tiles=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
    for j=1:2
        data=cases{j};
        label=sprintf('%s %g×%g Pitch',labels(j),data.cfg.periods);
        ax=nexttile(tiles);
        imagesc(ax,coord,coord,10*log10(max(data.spectralEnergy(:,:,2)/greenReference,1e-8)));
        colormap(ax,gray(256)); clim(ax,[-60,0]); colorbar(ax);
        title(ax,label+" | 550 nm, dB"); formatAngleAxes(ax,half);
        ax=nexttile(tiles);
        image(ax,coord,coord,data.rgb);
        title(ax,label+" | RGB代理"); formatAngleAxes(ax,half);
        ax=nexttile(tiles);
        image(ax,coord,coord,data.rgb);
        title(ax,label+" | 中心放大"); formatAngleAxes(ax,zoomHalf);
    end
    title(tiles,'真实Anode/BM · 平面波 · 620/550/450 nm · 共同曝光及对数显示');
    set(findall(fig,'Type','axes'),'FontName','Microsoft YaHei','FontSize',10);
    exportgraphics(fig,fullfile(outputDir,'plane_wave_comparison.png'),'Resolution',160);

    fig2=figure('Color','w','Position',[80,80,1100,500],'Name','平面波反射区域与角剖面');
    tiles=tiledlayout(fig2,1,2,'TileSpacing','compact','Padding','compact');
    ax=nexttile(tiles);
    n=size(report.local.effectiveReflection,1);
    dx=report.local.cfg.pitch/report.local.cfg.samplesPerPeriod;
    xy=((0:n-1)-floor(n/2))*dx*1e6;
    imagesc(ax,xy,xy,abs(report.local.effectiveReflection).^2);
    axis(ax,'image'); set(ax,'YDir','normal'); colormap(ax,gray(256)); clim(ax,[0,1]);
    xlabel(ax,'x / μm'); ylabel(ax,'y / μm'); title(ax,'原图完整单元：有效强度反射率');
    ax=nexttile(tiles); hold(ax,'on');
    mid=(numel(coord)+1)/2;
    for j=1:2
        curve=cases{j}.spectralEnergy(mid,:,2);
        plot(ax,coord,10*log10(max(curve/greenReference,1e-8)), ...
            'DisplayName',labels(j),'LineWidth',1);
    end
    xlabel(ax,'θx / °'); ylabel(ax,'550 nm / dB，共同参考');
    xlim(ax,[-zoomHalf,zoomHalf]); ylim(ax,[-60,0]); grid(ax,'on'); legend(ax,'Location','best');
    title(ax,'θy=0 剖面：同一角度像素的积分能量');
    set(findall(fig2,'Type','axes'),'FontName','Microsoft YaHei','FontSize',10);
    exportgraphics(fig2,fullfile(outputDir,'geometry_and_profile.png'),'Resolution',160);
end

function formatAngleAxes(ax,half)
    axis(ax,'image'); set(ax,'YDir','normal');
    xlim(ax,[-half,half]); ylim(ax,[-half,half]);
    xlabel(ax,'θx / °'); ylabel(ax,'θy / °');
end
