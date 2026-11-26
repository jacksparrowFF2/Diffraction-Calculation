function plotImagingComparison(optics,outputDir)
%plotImagingComparison 对比同一入瞳的无遮挡与含 Mask 响应。
    green = 2;
    x = optics.coordinate*1e6;
    clearPsf = optics.quality.clearPsf(:,:,green);
    maskPsf = optics.psf(:,:,green);
    referencePeak = max(clearPsf,[],'all');
    clearDb = 10*log10(max(clearPsf/referencePeak,1e-7));
    maskDb = 10*log10(max(maskPsf/referencePeak,1e-7));
    selected = abs(x)<=25+8*eps(25);
    if nnz(selected)<3, selected(:) = true; end
    detail = x(selected);

    fig = figure('Visible','off','Color','white','Position',[100,100,920,1040]);
    closeFigure = onCleanup(@() close(fig));
    layout = tiledlayout(fig,3,2,'TileSpacing','compact','Padding','compact');
    drawPsf(nexttile(layout),x,clearDb,'无遮挡 · 仅入瞳');
    ax = nexttile(layout);
    drawPsf(ax,x,maskDb,'含 UDC Mask');
    scale = colorbar(ax);
    scale.Layout.Tile = 'east';
    scale.Label.String = '相对无遮挡峰值 (dB)';
    drawPsf(nexttile(layout),detail,clearDb(selected,selected),'无遮挡 · 中心放大');
    drawPsf(nexttile(layout),detail,maskDb(selected,selected),'含 UDC Mask · 中心放大');

    profiles = optics.quality.mtfProfiles;
    profiles = profiles(profiles.Channel==green,:);
    frequency = unique(profiles.CyclesPerMm,'stable');
    masked = min(reshape(profiles.MTF,numel(frequency),[]),[],2);
    clearMtf = min(reshape(profiles.ClearMTF,numel(frequency),[]),[],2);
    ax = nexttile(layout,[1,2]);
    plot(ax,frequency,clearMtf,'--','Color',[0.39,0.45,0.53],'LineWidth',1.6,'DisplayName','无遮挡');
    hold(ax,'on');
    plot(ax,frequency,masked,'Color',[0.13,0.34,0.58],'LineWidth',1.8,'DisplayName','含 UDC Mask');
    title(ax,'G 通道 MTF · 逐频率取四方向最小值');
    xlabel(ax,'空间频率 (cycles/mm)');
    ylabel(ax,'MTF');
    ylim(ax,[0,1.02]);
    xlim(ax,[0,frequency(end)]);
    grid(ax,'on');
    legend(ax,'Location','best');
    title(layout,'同一入瞳下的 G 通道归一化 PSF 与 MTF');
    exportgraphics(layout,fullfile(outputDir,'udc_imaging_comparison.png'),'Resolution',160);
end

function drawPsf(ax,x,db,heading)
    imagesc(ax,x,x,db);
    axis(ax,'image');
    set(ax,'YDir','normal');
    colormap(ax,turbo);
    clim(ax,[-70,0]);
    title(ax,heading);
    xlabel(ax,'传感器 x (μm)');
    ylabel(ax,'传感器 y (μm)');
end
