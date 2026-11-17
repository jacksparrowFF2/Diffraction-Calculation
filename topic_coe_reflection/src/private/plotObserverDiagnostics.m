function [regionFigure,stageFigure]=plotObserverDiagnostics(views,count,rgbOptions)
%plotObserverDiagnostics 几何贡献区域，以及点源/有限源/入瞳三阶段。
    peak=0; extent=0;
    for j=1:numel(views)
        peak=max(peak,max(views{j}.participation.rayPairDensity,[],'all'));
        extent=max(extent,max(abs(views{j}.participation.coordinate)));
    end
    regionFigure=figure('Name','镜面几何贡献区域','Color','w','Position',[40,70,1500,760]);
    layout=tiledlayout(2,count,'TileSpacing','compact','Padding','compact');
    title(layout,sprintf('镜面零级几何贡献区域 · 模型%d×%d Pitch · 上排移动光源，下排移动观察者', ...
        views{1}.config.periods(1),views{1}.config.periods(2)));
    subtitle(layout,'灰度：光线对密度，共同尺度；青圈：理论支撑范围；绿框：实际计算面；不裁切衍射场');
    for j=1:numel(views)
        nexttile; view=views{j}; region=view.participation;
        imagesc(region.coordinate*1e3,region.coordinate*1e3,region.rayPairDensity/max(peak,realmin));
        axis image xy; colormap(gca,gray(256)); clim([0,1]); hold on;
        t=linspace(0,2*pi,361); r=region.diameter/2;
        plot(r*cos(t)*1e3,r*sin(t)*1e3,'c-','LineWidth',1);
        bounds=region.inputBounds*1e3;
        rectangle('Position',[bounds(1),bounds(3),diff(bounds(1:2)),diff(bounds(3:4))], ...
            'EdgeColor',[0.2,0.9,0.2],'LineStyle','--');
        xlim([-extent,extent]*1e3); ylim([-extent,extent]*1e3);
        xlabel('屏幕 x / mm'); ylabel('屏幕 y / mm');
        title(sprintf('a=%.0f, b=%.0f cm\n几何D=%.3f mm，跨径%.1f Pitch', ...
            100*view.sourceDistance,100*view.observerDistance,region.diameter*1e3,region.nominalPeriods));
    end
    distances=cellfun(@(v) (v.sourceDistance-0.25)^2+(v.observerDistance-0.25)^2,views);
    [~,chosen]=min(distances); view=views{chosen};
    weights=view.config.channelResponse.*view.config.spectralWeights(:);
    arrays={view.pointSpectralEnergy,view.sourceAveragedSpectralEnergy,view.spectralEnergy};
    images=cell(1,3); n=numel(view.coordinate);
    for j=1:3
        images{j}=reshape(reshape(arrays{j},[],size(weights,1))*weights,n,n,3)/view.binSolidAngle;
    end
    pooled=vertcat(images{:}); [rgb,~,~]=plotComposeReflectionRgb(pooled,pooled,rgbOptions,false);
    stageFigure=figure('Name','有限光源和入瞳作用分解','Color','w','Position',[60,80,1200,760]);
    layout=tiledlayout(2,3,'TileSpacing','compact');
    title(layout,sprintf('计算%d×%d Pitch；a=%.0f cm, b=%.0f cm；源%.1f mm，入瞳%.1f mm；三图共同曝光', ...
        view.config.periods(1),view.config.periods(2),100*view.sourceDistance, ...
        100*view.observerDistance,view.sourceDiameter*1e3,view.pupilDiameter*1e3));
    labels=["点源参考","均匀圆盘源强度平均","再加理想入瞳强度PSF（可选）"];
    coordinate=rad2deg(view.coordinate);
    for j=1:3
        nexttile; rows=(j-1)*n+(1:n); imagesc(coordinate,coordinate,rgb(rows,:,:));
        axis image xy; xlabel('观察角 / °'); ylabel('观察角 / °'); title(labels(j));
    end
    nexttile([1,3]); [~,k]=min(abs(view.config.wavelengths-550e-9));
    centre=(n+1)/2; reference=max(arrays{1}(:,:,k),[],'all');
    for j=1:3
        plot(coordinate,10*log10(max(arrays{j}(centre,:,k)/max(reference,realmin),1e-10)),'LineWidth',1.2); hold on;
    end
    grid on; xlim(coordinate([1,end])); ylim([-70,1]); xlabel('观察角 theta_x / °'); ylabel('共同参考 / dB');
    legend(labels,'Location','best'); title(sprintf('%.0f nm 中心截线；不对各阶段单独归一化',view.config.wavelengths(k)*1e9));
end
