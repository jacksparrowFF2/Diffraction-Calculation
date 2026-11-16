function report=benchFourPatterns(cfg,stack,inputRoot,outputRoot,rgbOptions)
%benchFourPatterns 用户指定四图样：各局部图样及其有限重复阵列，共同曝光RGB。
% 图样由Input/case01~case04读取；正式计算不生成、不更改轮廓。
    names=["case01","case02","case03","case04"];
    labels=["单圆","等面积椭圆","四圆阵列","四椭圆（不同取向）"];
    report.protocol="理想二值反射图样；平面波；620/550/450nm RGB代理；所有图共同曝光；无相机模型";
    report.caseNames=names; report.labels=labels;
    report.patterns=cell(1,4); report.local=cell(1,4); report.array=cell(1,4);
    summary=table(); started=tic;
    for j=1:4
        inputDir=fullfile(inputRoot,names(j));
        settings=struct('cfg',cfg,'stack',stack,'spectrumMode',"RGB", ...
            'rgbOptions',rgbOptions,'previewWavelengthNm',[620,550,450], ...
            'spotOptions',struct('wavelengthNm',[620,550,450]), ...
            'observation',struct('centerDeg',[0,0],'halfAngleDeg',2.5, ...
                'innerHalfAngleDeg',0,'curveSamples',301), ...
            'displayStyle',struct('showAnnotations',false), ...
            'saveResults',true,'showProgress',true,'caseName',"coe_reflection");
        [settings,record]=coreLoadProjectConfig(fullfile(inputDir,'config.txt'),settings);
        current=settings.cfg; currentStack=settings.stack;
        assert(current.source.type=="plane" && all(current.heightMap(:)==0) ...
            && isempty(current.complexCoefficient) && isinf(current.pupilRadius) ...
            && isequal(current.channelResponse,eye(3)) && settings.spectrumMode=="RGB", ...
            'COEReflection:ComparisonInput','四图样对比要求平面波、零额外高度、图片及RGB代理。');
        assert(numel(current.wavelengths)==3 && max(abs(current.wavelengths-[620,550,450]*1e-9))<1e-15 ...
            && isequal(currentStack.layerNames,["Anode","BM"]), ...
            'COEReflection:ComparisonInput','要求Anode/BM和620/550/450nm。');
        [localCfg,native]=benchInputCell(current,currentStack,inputDir);
        assert(native==current.samplesPerPeriod && any(current.periods>localCfg.periods), ...
            'COEReflection:ComparisonSampling','保留原生采样，目标阵列须大于局部单元。');
        variants={localCfg,current};
        for state=1:2
            fprintf('四图样 %s %d/2：%g×%g Pitch\n',names(j),state,variants{state}.periods);
            value=runCOEReflection(variants{state},currentStack,inputDir);
            assert(value.propagationRuns==1 && max(value.metrics.ParsevalAbsoluteError)<1e-10, ...
                'COEReflection:ComparisonEnergy','传播次数或能量检查失败。');
            data=struct('cfg',variants{state},'coordinate',value.coordinate, ...
                'spectralEnergy',value.spectralEnergy,'metrics',value.metrics, ...
                'incidentArea',value.incidentArea,'incidentPower',value.incidentPower, ...
                'illuminationAperture',value.illuminationAperture, ...
                'repeatCounts',value.layerGeometry.repeatCounts,'projectConfig',record);
            if state==1
                model=value.config.coeStack;
                reflection=reshape(model.coefficientLookup(double(model.regionMap)+1,2),size(model.regionMap));
                report.patterns{j}=abs(reflection).^2;
                report.local{j}=data;
            else
                report.array{j}=data;
                assert(data.illuminationAperture.coverageFraction<1 ...
                    || report.local{j}.illuminationAperture.coverageFraction<1 ...
                    || max(abs(data.metrics.Throughput-report.local{j}.metrics.Throughput))<1e-12, ...
                    'COEReflection:ComparisonTiling','完整平铺前后总反射率不同。');
            end
            row=data.metrics; row.Case=repmat(names(j),height(row),1);
            row.IsArray=repmat(state==2,height(row),1);
            row.PeriodsX=repmat(data.cfg.periods(1),height(row),1);
            row.PeriodsY=repmat(data.cfg.periods(2),height(row),1);
            summary=[summary;row]; %#ok<AGROW>
            clear value model reflection
        end
    end
    % 统一角坐标和波长权重；不分别按case调整曝光或色彩。
    for j=1:4
        for data={report.local{j},report.array{j}}
            assert(isequal(data{1}.coordinate,report.local{1}.coordinate) ...
                && isequal(data{1}.cfg.spectralWeights,report.local{1}.cfg.spectralWeights), ...
                'COEReflection:ComparisonGrid','各case须有共同角网格和光谱权重。');
        end
    end
    weights=reshape(report.local{1}.cfg.spectralWeights,1,1,3);
    localSignal=cellfun(@(x) x.spectralEnergy.*weights,report.local,'UniformOutput',false);
    arraySignal=cellfun(@(x) x.spectralEnergy.*weights,report.array,'UniformOutput',false);
    [localRGB,arrayRGB,report.mapping]=plotComposeReflectionRgb( ...
        vertcat(localSignal{:}),vertcat(arraySignal{:}),rgbOptions,true);
    clear localSignal arraySignal
    count=numel(report.local{1}.coordinate);
    localImages=cell(1,4); arrayImages=cell(1,4);
    for j=1:4
        rows=(j-1)*count+(1:count);
        localImages{j}=localRGB(rows,:,:); arrayImages{j}=arrayRGB(rows,:,:);
        folder=fullfile(outputRoot,names(j)); if ~isfolder(folder), mkdir(folder); end
        imwrite(flipud(localImages{j}),fullfile(folder,'local_rgb.png'));
        imwrite(flipud(arrayImages{j}),fullfile(folder,'array_rgb.png'));
        writetable(summary(summary.Case==names(j),:),fullfile(folder,'pattern_metrics.csv'));
    end
    report.summary=summary; report.elapsedSeconds=toc(started);
    outputDir=fullfile(outputRoot,'FourPatterns'); if ~isfolder(outputDir), mkdir(outputDir); end
    save(fullfile(outputDir,'four_patterns_results.mat'),'report','-v7.3');
    writetable(summary,fullfile(outputDir,'four_patterns_summary.csv'));
    renderPatterns(report,localImages,arrayImages,outputDir);
    fprintf('四图样对比完成，%.1f s：%s\n',toc(started),outputDir);
end

function renderPatterns(report,localImages,arrayImages,outputDir)
    coord=rad2deg(report.local{1}.coordinate);
    reference=0;
    for j=1:4
        reference=max(reference,max(report.local{j}.spectralEnergy(:,:,2),[],'all'));
    end
    fig=figure('Color','w','Position',[50,30,1650,1150],'Name','case01~case04平面波RGB对比');
    tiles=tiledlayout(fig,4,4,'TileSpacing','compact','Padding','compact');
    for row=1:4
        for j=1:4
            ax=nexttile(tiles);
            if row==1
                pattern=report.patterns{j};
                xy=((0:size(pattern,1)-1)-size(pattern,1)/2) ...
                    *report.local{j}.cfg.pitch/report.local{j}.cfg.samplesPerPeriod*1e6;
                imagesc(ax,xy,xy,pattern); colormap(ax,gray(256)); clim(ax,[0,1]);
                xlabel(ax,'x / μm'); ylabel(ax,'y / μm');
                title(ax,report.caseNames(j)+" "+report.labels(j));
            elseif row==2
                intensity=report.local{j}.spectralEnergy(:,:,2);
                imagesc(ax,coord,coord,10*log10(max(intensity/reference,1e-8)));
                colormap(ax,gray(256)); clim(ax,[-60,0]);
                title(ax,'局部图样 550 nm / dB，共同参考');
            elseif row==3
                image(ax,coord,coord,localImages{j}); title(ax,'局部图样 RGB代理');
            else
                image(ax,coord,coord,arrayImages{j}); title(ax,'重复阵列 8×8 Pitch RGB代理');
            end
            axis(ax,'image'); set(ax,'YDir','normal');
            if row>1
                half=report.local{j}.cfg.angleHalfRangeDeg;
                xlim(ax,[-half,half]); ylim(ax,[-half,half]);
                xlabel(ax,'θx / °'); ylabel(ax,'θy / °');
            end
        end
    end
    title(tiles,'四种理想反射图样 · 平面波 · 620/550/450 nm · 全部RGB共用曝光');
    set(findall(fig,'Type','axes'),'FontName','Microsoft YaHei','FontSize',9);
    exportgraphics(fig,fullfile(outputDir,'four_patterns_comparison.png'),'Resolution',160);
    fig2=figure('Color','w','Position',[70,60,1200,620],'Name','四图样局部RGB');
    tiles=tiledlayout(fig2,1,2,'TileSpacing','compact','Padding','compact');
    for j=3:4
        ax=nexttile(tiles); image(ax,coord,coord,localImages{j});
        axis(ax,'image'); set(ax,'YDir','normal');
        half=report.local{j}.cfg.angleHalfRangeDeg;
        xlim(ax,[-half,half]); ylim(ax,[-half,half]);
        xlabel(ax,'θx / °'); ylabel(ax,'θy / °');
        title(ax,report.caseNames(j)+" "+report.labels(j));
    end
    title(tiles,'四圆与四椭圆局部图样：共同曝光RGB代理');
    set(findall(fig2,'Type','axes'),'FontName','Microsoft YaHei','FontSize',11);
    exportgraphics(fig2,fullfile(outputDir,'case03_case04_local_rgb.png'),'Resolution',160);
end
