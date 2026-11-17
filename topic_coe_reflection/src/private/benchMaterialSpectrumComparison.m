function report = benchMaterialSpectrumComparison(defaultSettings, moduleDir, outputDir)
%benchMaterialSpectrumComparison 三个实际项目的全谱正式传播与统一曝光比较。
% 逐项目计算/释放大数组；原始能量为入射功率归一化的double，不保存可推导的重复数组。
    if nargin < 3
        stamp = string(datetime('now','TimeZone','Asia/Shanghai','Format','yyyyMMdd_HHmmss'));
        outputDir = fullfile(moduleDir,'Output',"MaterialSpectra_"+stamp);
    end
    if ~isfolder(outputDir), mkdir(outputDir); end
    diary(fullfile(outputDir,'simulation.log'));
    logCleanup = onCleanup(@() diary('off')); %#ok<NASGU>
    names = ["F192","F192Base","COEB3"];
    report = struct('outputDir',string(outputDir),'projects',names,'referencePeak',0);
    summary = table(); spectralSummary = table(); roiSummary = table();
    for j = 1:numel(names)
        name = names(j); inputDir = fullfile(moduleDir,'Input',name);
        projectDir = fullfile(outputDir,name);
        if ~isfolder(projectDir), mkdir(projectDir); end
        [settings,projectConfig] = coreLoadProjectConfig(fullfile(inputDir,'config.txt'),defaultSettings);
        cfg = settings.cfg; stack = settings.stack;
        assert(settings.spectrumMode=="D65" && settings.wavelengthStepNm==10 ...
            && isequal(settings.wavelengthRangeNm,[380,780]) ...
            && isequal(cfg.periods,[4,4]) && cfg.samplesPerPeriod==64 ...
            && cfg.outputSize==1801 && cfg.angleHalfRangeDeg==5 ...
            && cfg.propagationDistance==0.25 && cfg.incidentIrradiance==1 ...
            && cfg.source.type=="plane" && cfg.source.integrationMethod=="direct" ...
            && cfg.source.apertureShape=="none" && isinf(cfg.pupilRadius) ...
            && ~cfg.observer.enabled && all(cfg.heightMap(:)==0), ...
            'COEReflection:MaterialProtocol','项目参数不符合正式全谱比较协议。');
        spectrum = coreLoadSpectrum(380:10:780,fullfile(moduleDir,'Input','Colorimetry'));
        cfg.wavelengths = spectrum.wavelengthNm.'/1e9;
        cfg.spectralWeights = spectrum.energyWeights.';
        cfg.channelResponse = ones(41,1); cfg.observation = settings.observation;
        for layerName = ["BM","CF_R","CF_G","CF_B"]
            layer = stack.layers.(layerName);
            assert(layer.method=="nk" && layer.nkFile=="Materials/"+layerName+"_nk.csv", ...
                'COEReflection:MaterialPath','未启用预期材料文件。');
            nk = readmatrix(fullfile(inputDir,layer.nkFile));
            assert(size(nk,2)==3 && all(isfinite(nk),'all') && all(diff(nk(:,1))>0) ...
                && all(nk(:,2)>0) && all(nk(:,3)>=0) && nk(1,1)<=380 && nk(end,1)>=780, ...
                'COEReflection:MaterialData','材料文件无效或未覆盖计算范围。');
        end
        started = tic;
        fprintf('\n===== %s: 正式41波长 × 1801²，4×4 Pitch =====\n',name);
        value = runCOEReflection(cfg,stack,inputDir);
        cfg = value.config; geometry = value.layerGeometry;
        assert(value.propagationRuns==1 && isequal(size(value.spectralEnergy),[1801,1801,41]) ...
            && all(isfinite(value.spectralEnergy),'all') && all(value.spectralEnergy>=0,'all') ...
            && max(value.metrics.ParsevalAbsoluteError)<1e-10 ...
            && all(value.metrics.CapturedThroughput<=value.metrics.Throughput+1e-8), ...
            'COEReflection:MaterialNumerics','完整传播或能量检查失败。');
        verifyLoadedMaterials(cfg,stack,inputDir);
        phaseStep = warnThickFilm(cfg);
        angular = coreComputeAngularMetrics(value,settings.observation,spectrum);
        spots = coreComputeDiffractionSpots(value,value,spectrum,settings.spotOptions,false);
        [~,color,mapping] = coreSynthesizeColor(value,value,spectrum,settings.colorOptions,false);
        report.referencePeak = max(report.referencePeak,mapping.referencePeak);
        windowXYZ = squeeze(sum(sum(color.XYZ,1),2)).';
        windowXY = windowXYZ(1:2)/sum(windowXYZ);
        [xx,yy] = meshgrid(rad2deg(value.coordinate));
        cone = hypot(xx,yy)<=2.5;
        viewXYZ = sum(reshape(color.XYZ,[],3).*cone(:),1);
        viewXY = viewXYZ(1:2)/sum(viewXYZ);
        mid = (cfg.outputSize+1)/2;
        [~,green] = min(abs(cfg.wavelengths-550e-9));
        greenDensity = value.spectralEnergy(mid,:,green)/(value.coordinate(2)-value.coordinate(1))^2;
        photopicDensity = color.XYZDensity(mid,:,2);
        profiles = table(rad2deg(value.coordinate(:)),greenDensity(:),photopicDensity(:), ...
            'VariableNames',{'ThetaXDeg','GreenDensityPerSr','PhotopicDensityPerSr'});
        greenFwhm = centralFwhm(profiles.ThetaXDeg,greenDensity(:));
        photoFwhm = centralFwhm(profiles.ThetaXDeg,photopicDensity(:));
        elapsedSeconds = toc(started);
        row = table(name,cfg.pitch*1e6,value.incidentArea,value.incidentPower, ...
            100*angular.radiantTotalReflectance,100*angular.photopicTotalReflectance, ...
            100*angular.radiantCapturedReflectance,100*angular.photopicCapturedReflectance, ...
            100*angular.radiantCapturedReflectance/angular.radiantTotalReflectance, ...
            100*angular.photopicCapturedReflectance/angular.photopicTotalReflectance, ...
            100*angular.radiantViewReflectance,100*angular.photopicViewReflectance, ...
            windowXY(1),windowXY(2),viewXY(1),viewXY(2), ...
            rad2deg(550e-9/cfg.pitch),cfg.propagationDistance*550e-9/cfg.pitch*1e3, ...
            greenFwhm,photoFwhm,value.fraunhoferNumber,phaseStep, ...
            max(value.metrics.ParsevalAbsoluteError),elapsedSeconds, ...
            'VariableNames',{'Project','PitchUm','IncidentAreaM2','IncidentPowerW', ...
            'RadiantTotalPercent','PhotopicTotalPercent','RadiantWindowPercent','PhotopicWindowPercent', ...
            'RadiantCaptureOfReflectionPercent','PhotopicCaptureOfReflectionPercent', ...
            'RadiantCone2p5Percent','PhotopicCone2p5Percent','WindowCIE_x','WindowCIE_y', ...
            'ConeCIE_x','ConeCIE_y','Nominal550PitchOrderDeg','Nominal550PitchOrderMm', ...
            'GreenCentralFwhmDeg','PhotopicCentralFwhmDeg','FraunhoferNumber', ...
            'MaxRoundTripPhaseStepRad','MaxParsevalError','SimulationSeconds'});
        zeroAbsolute = value.metrics.Throughput.*value.metrics.ZeroOrderFraction;
        row.RadiantZeroPercent = 100*dot(zeroAbsolute,spectrum.energyWeights);
        row.PhotopicZeroPercent = 100*dot(zeroAbsolute,spectrum.photopicWeights);
        row.RadiantZeroOfReflectionPercent = 100*row.RadiantZeroPercent/row.RadiantTotalPercent;
        row.PhotopicZeroOfReflectionPercent = 100*row.PhotopicZeroPercent/row.PhotopicTotalPercent;
        row.RadiantOutside2p5Percent = row.RadiantTotalPercent-row.RadiantCone2p5Percent;
        row.PhotopicOutside2p5Percent = row.PhotopicTotalPercent-row.PhotopicCone2p5Percent;
        summary = [summary;row]; %#ok<AGROW>
        metrics = value.metrics; metrics.Project = repmat(name,height(metrics),1);
        spectralSummary = [spectralSummary;metrics]; %#ok<AGROW>
        for radius = [0.15,0.30,0.5,1,2.5,5]
            roi = hypot(xx,yy)<=radius;
            perLambda = sum(reshape(value.spectralEnergy,[],41).*roi(:),1).';
            roiSummary = [roiSummary;table(name,radius,100*dot(perLambda,spectrum.energyWeights), ...
                100*dot(perLambda,spectrum.photopicWeights),'VariableNames', ...
                {'Project','RadiusDeg','RadiantIncidentPercent','PhotopicIncidentPercent'})]; %#ok<AGROW>
        end
        physical = struct('spectralEnergy',value.spectralEnergy,'coordinate',value.coordinate, ...
            'receiverCoordinate',value.receiverCoordinate,'receiverPixelArea',value.receiverPixelArea, ...
            'incidentArea',value.incidentArea,'incidentPower',value.incidentPower, ...
            'metrics',value.metrics,'propagationMethod',value.propagationMethod, ...
            'receiverSamplesPerBin',value.receiverSamplesPerBin,'fraunhoferNumber',value.fraunhoferNumber);
        storage = "spectralEnergy: per-bin reflected power / incident power, double; " + ...
            "spectralIrradiance = spectralEnergy * incidentPower / receiverPixelArea (W/m²); " + ...
            "no duplicate baseline/clear/irradiance arrays; XYZ uses saved spectrum weights";
        save(fullfile(projectDir,'physical.mat'),'physical','cfg','stack','geometry', ...
            'spectrum','settings','projectConfig','angular','spots','storage','-v7.3');
        writelines(projectConfig.text,fullfile(projectDir,'config_snapshot.txt'),'Encoding','UTF-8');
        writetable(cfg.coeStack.materialTable,fullfile(projectDir,'stack_materials.csv'));
        writetable(cfg.coeStack.regionTable,fullfile(projectDir,'stack_regions.csv'));
        writetable(metrics,fullfile(projectDir,'spectral_metrics.csv'));
        writetable(spots.metrics,fullfile(projectDir,'diffraction_spots.csv'));
        writetable(profiles,fullfile(projectDir,'horizontal_profiles.csv'));
        cumulative = table(angular.radiusDeg,100*angular.radiantCumulativeReflectance, ...
            100*angular.photopicCumulativeReflectance,'VariableNames', ...
            {'RadiusDeg','RadiantIncidentPercent','PhotopicIncidentPercent'});
        writetable(cumulative,fullfile(projectDir,'cumulative_energy.csv'));
        plotGeometry(geometry,stack,cfg,name,projectDir);
        fprintf('%s 完成：总反射 %.6f%%，明视觉 %.6f%%；计算 %.1fs。\n', ...
            name,row.RadiantTotalPercent,row.PhotopicTotalPercent,elapsedSeconds);
        writetable(summary,fullfile(outputDir,'comparison_summary.csv'));
        clear value physical color angular spots xx yy cone roi cfg geometry
    end
    report.fixedGain = 0.9/report.referencePeak;
    report.summary = summary;
    report.protocol = "same incident irradiance 1 W/m²; spectralEnergy/incident power; " + ...
        "D65 380:10:780 nm trapezoid; shared gain = 0.9 / largest linear RGB density peak across all projects; " + ...
        "same EV=0, toneStrength=60, luminanceGamma=0.72; no centre suppression";
    writetable(spectralSummary,fullfile(outputDir,'spectral_metrics_all.csv'));
    writetable(roiSummary,fullfile(outputDir,'fixed_roi_comparison.csv'));
    writetable(table(spectrum.wavelengthNm,spectrum.incidentSPD,spectrum.quadratureNm, ...
        spectrum.energyWeights,spectrum.photopicWeights,spectrum.cmfXYZ, ...
        'VariableNames',{'WavelengthNm','D65SPD','QuadratureNm','EnergyWeight','PhotopicWeight','CMF_XYZ'}), ...
        fullfile(outputDir,'spectrum_weights.csv'));
    fprintf('\n跨项目共同曝光 fixedGain=%.12g；统一映射开始。\n',report.fixedGain);
    renderAll(report);
    save(fullfile(outputDir,'comparison_report.mat'),'report','spectralSummary','roiSummary');
    fprintf('全部完成：%s\n',outputDir);
end

function verifyLoadedMaterials(cfg,stack,inputDir)
    materials = cfg.coeStack.materialTable;
    for name = ["BM","CF_R","CF_G","CF_B"]
        data = readmatrix(fullfile(inputDir,stack.layers.(name).nkFile));
        chosen = materials.Material==name;
        n = interp1(data(:,1),data(:,2),materials.WavelengthNm(chosen));
        k = interp1(data(:,1),data(:,3),materials.WavelengthNm(chosen));
        assert(nnz(chosen)==41 && max(abs(materials.n(chosen)-n))<1e-12 ...
            && max(abs(materials.k(chosen)-k))<1e-12 ...
            && all(materials.ThicknessM(chosen)==stack.layers.(name).thickness), ...
            'COEReflection:LoadedMaterial','实际膜堆没有按原膜厚加载材料CSV。');
    end
end

function phaseStep = warnThickFilm(cfg)
    materials = cfg.coeStack.materialTable; opticalDepth = 0;
    for level = unique(materials.Level).'
        selected = materials.Level==level & materials.Method=="nk";
        if any(selected)
            opticalDepth = opticalDepth+max(materials.n(selected).*materials.ThicknessM(selected));
        end
    end
    phaseStep = 4*pi*opticalDepth*max(abs(diff(1./cfg.wavelengths)));
    if phaseStep>pi/2
        warning('COEReflection:ThickFilmSpectrum', ...
            '厚膜光谱干涉可能欠采样（最大往返相位步进约 %.2f rad）；10nm结果尚无光谱收敛保证，建议步长约 %.3g nm或更小并检查收敛。', ...
            phaseStep,min(cfg.wavelengths)^2/(8*opticalDepth)*1e9);
    end
end

function width = centralFwhm(x,y)
    mid = (numel(x)+1)/2; width = NaN;
    if y(mid)<=0 || y(mid)<max(y(mid-1:mid+1)), return; end
    half = y(mid)/2;
    left = find(y(1:mid)<half,1,'last');
    right = mid-1+find(y(mid:end)<half,1,'first');
    if isempty(left)||isempty(right), return; end
    a = interp1(y(left:left+1),x(left:left+1),half);
    b = interp1(y(right-1:right),x(right-1:right),half);
    width = b-a;
end

function plotGeometry(geometry,stack,cfg,name,projectDir)
    fig = figure('Visible','off','Color','w','Position',[40,40,1250,700]);
    tiles = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
    xy = ((0:255)-128)*cfg.pitch/cfg.samplesPerPeriod*1e6;
    for layerName = stack.layerNames
        ax = nexttile(tiles); imagesc(ax,xy,xy,geometry.layerMasks.(layerName));
        axis(ax,'image'); set(ax,'YDir','normal'); clim(ax,[0,1]); colormap(ax,gray(256));
        title(ax,layerName,'Interpreter','none'); xlabel(ax,'x / μm'); ylabel(ax,'y / μm');
    end
    title(tiles,sprintf('%s | 4×4 Pitch | 64 samples/Pitch | native 250 samples/Pitch',name));
    exportgraphics(fig,fullfile(projectDir,'sampled_geometry.png'),'Resolution',130); close(fig);
end

function renderAll(report)
    enhanced = cell(1,3); linear = cell(1,3);
    for j = 1:3
        name = report.projects(j); projectDir = fullfile(report.outputDir,name);
        data = load(fullfile(projectDir,'physical.mat'),'physical','cfg','spectrum','settings');
        value = struct('config',data.cfg,'coordinate',data.physical.coordinate, ...
            'spectralEnergy',data.physical.spectralEnergy);
        options = data.settings.colorOptions;
        assert(options.toneStrength==60 && options.luminanceGamma==0.72 && options.exposureEV==0, ...
            'COEReflection:MaterialDisplay','增强显示参数与协议不符。');
        options.fixedGain = report.fixedGain;
        [~,color,mapping] = coreSynthesizeColor(value,value,data.spectrum,options,false);
        enhanced{j} = color.displaySRGB; linear{j} = color.sRGB;
        xyzDensity = color.XYZDensity; linearSRGBDensity = color.linearSRGBDensity;
        save(fullfile(projectDir,'color_physical.mat'),'xyzDensity','linearSRGBDensity','mapping','-v7.3');
        imwrite(flipud(enhanced{j}),fullfile(projectDir,'rgb_enhanced_shared.png'));
        imwrite(flipud(linear{j}),fullfile(projectDir,'rgb_linear_shared.png'));
        coord = rad2deg(value.coordinate); mm = value.coordinate*250;
        for mode = ["enhanced","linear"]
            rgb = enhanced{j}; if mode=="linear", rgb=linear{j}; end
            exportRgb(rgb,coord,name,mode,"angle",projectDir,report.fixedGain);
            exportRgb(rgb,mm,name,mode,"receiver_mm",projectDir,report.fixedGain);
        end
        % 保留既有毫米坐标工具的距离/功率表；正式RGB图以上方完整标题与统一曝光为准。
        receiver = data.physical; receiver.config = data.cfg;
        receiver.spectralIrradiance = receiver.spectralEnergy*receiver.incidentPower/receiver.receiverPixelArea;
        oldVisible = get(groot,'DefaultFigureVisible');
        visibilityCleanup = onCleanup(@() set(groot,'DefaultFigureVisible',oldVisible));
        set(groot,'DefaultFigureVisible','off');
        [f1,f2,view] = plotReceiverDistanceViews(receiver,enhanced{j},0.25);
        view.displayConvention = "cross-project shared exposure and mapping; no per-project normalization";
        writetable(view.summary,fullfile(projectDir,'receiver_distance_summary.csv'));
        close(f1); close(f2); clear visibilityCleanup
        clear data value color receiver xyzDensity linearSRGBDensity
    end
    for mode = ["enhanced","linear"]
        images = enhanced; if mode=="linear", images=linear; end
        for units = ["angle","receiver_mm"]
            coordinate = coord; if units=="receiver_mm", coordinate=mm; end
            fig = figure('Visible','off','Color','w','Position',[20,20,1800,690]);
            tiles = tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
            for j = 1:3
                ax = nexttile(tiles); image(ax,coordinate,coordinate,images{j});
                formatRgbAxes(ax,units); title(ax,report.projects(j),'Interpreter','none');
            end
            title(tiles,'F192 / F192Base / COEB3 | 4×4 Pitch | z=0.25 m | ±5° | D65 380:10:780 nm');
            subtitle(tiles,displayLabel(mode,report.fixedGain));
            exportgraphics(fig,fullfile(report.outputDir,"comparison_"+units+"_"+mode+".png"),'Resolution',150);
            close(fig);
        end
    end
end

function exportRgb(rgb,coordinate,name,mode,units,projectDir,gain)
    fig = figure('Visible','off','Color','w','Position',[40,40,1040,1000]);
    ax = axes(fig); image(ax,coordinate,coordinate,rgb); formatRgbAxes(ax,units);
    title(ax,sprintf('%s | 4×4 Pitch | z=0.25 m | ±5° | D65 380:10:780 nm',name),'Interpreter','none');
    subtitle(ax,displayLabel(mode,gain));
    exportgraphics(fig,fullfile(projectDir,"rgb_"+units+"_"+mode+"_shared.png"),'Resolution',150);
    close(fig);
end

function label = displayLabel(mode,gain)
    label = sprintf('Shared exposure gain=%.7g | linear exposure → sRGB',gain);
    if mode=="enhanced"
        label = sprintf('Shared exposure gain=%.7g | luminance log strength=60, gamma=0.72 → sRGB',gain);
    end
end

function formatRgbAxes(ax,units)
    axis(ax,'image'); set(ax,'YDir','normal','Color','k','FontSize',11);
    if units=="angle"
        xlabel(ax,'θx / °'); ylabel(ax,'θy / °');
    else
        xlabel(ax,'Receiver X / mm (X=zθx)'); ylabel(ax,'Receiver Y / mm (Y=zθy)');
    end
end
