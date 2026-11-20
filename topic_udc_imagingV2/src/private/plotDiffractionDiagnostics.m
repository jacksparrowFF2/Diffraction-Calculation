function plotDiffractionDiagnostics(optics,outputDir)
%plotDiffractionDiagnostics 导出有效入瞳、PSF、MTF 和同标尺点光源图。
    figureHandle = figure('Visible','off','Color','white','Position',[100,100,1440,390]);
    closeFigure = onCleanup(@() close(figureHandle));
    layout = tiledlayout(figureHandle,1,4,'TileSpacing','compact','Padding','compact');
    [~,wave] = min(abs(optics.config.wavelengths-550e-9));
    transmission = optics.transmissionLookup(:,wave);
    % 仅对预览抽样，避免把完整物面光学网格复制到绘图对象。
    stride = max(1,ceil(max(size(optics.pupil))/1500));
    cols = 1:stride:numel(optics.x);
    rows = 1:stride:numel(optics.y);
    region = optics.regionMap(rows,cols);
    amplitude = reshape(abs(transmission(double(region(:))+1)),size(region));
    ax = nexttile(layout);
    drawMap(ax,optics.x(cols)*1e6,optics.y(rows)*1e6,optics.pupil(rows,cols).*amplitude,'有效振幅入瞳','物面 (μm)');
    colormap(ax,gray); clim(ax,[0,1]); colorbar(ax);

    coordinate = optics.coordinate*1e6;
    ax = nexttile(layout);
    reference = max(max(optics.clearKernel,[],'all'),realmin);
    preview = min(log1p(1000*optics.kernel/reference)/log1p(1000),1);
    drawMap(ax,coordinate,coordinate,preview,'RGB 衍射强度（对数显示）','传感器 (μm)');

    ax = nexttile(layout);
    greenPsf = optics.psf(:,:,2);
    greenPeak = max(max(optics.quality.clearPsf(:,:,2),[],'all'),realmin);
    logPsf = 10*log10(max(greenPsf/greenPeak,1e-6));
    drawMap(ax,coordinate,coordinate,logPsf,'G 通道 PSF (dB)','传感器 (μm)');
    colormap(ax,turbo); clim(ax,[-60,0]); colorbar(ax);

    ax = nexttile(layout);
    frequency = optics.spatialFrequency;
    drawMap(ax,frequency,frequency,optics.mtf(:,:,2),'G 通道二维 MTF','空间频率 (cycles/mm)');
    colormap(ax,turbo); clim(ax,[0,1]); colorbar(ax);
    exportgraphics(layout,fullfile(outputDir,'udc_diffraction_diagnostics.png'),'Resolution',160);
    clear closeFigure

    figureHandle = figure('Visible','off','Color','white','Position',[100,100,900,460]);
    closeFigure = onCleanup(@() close(figureHandle));
    layout = tiledlayout(figureHandle,1,2,'TileSpacing','compact','Padding','compact');
    point = optics.quality.hdrPoint;
    reference = max([1,max(point.clear,[],'all'),max(point.withMask,[],'all')]);
    displayScale = log1p(20*reference);
    selected = abs(coordinate)<=optics.quality.options.hdrDisplayHalfWidthUm;
    displayCoordinate = coordinate(selected);
    drawMap(nexttile(layout),displayCoordinate,displayCoordinate,log1p(20*point.clear(selected,selected,:))/displayScale, ...
        '无遮挡','传感器 (μm)');
    drawMap(nexttile(layout),displayCoordinate,displayCoordinate,log1p(20*point.withMask(selected,selected,:))/displayScale, ...
        '含 UDC 掩膜','传感器 (μm)');
    title(layout,sprintf('HDR 点光源 · 输入 %.3g · 共同曝光与对数显示',optics.quality.options.hdrPointSignal));
    exportgraphics(layout,fullfile(outputDir,'udc_hdr_point.png'),'Resolution',160);
end

function drawMap(ax,x,y,data,label,unit)
    imagesc(ax,x,y,data);
    axis(ax,'image');
    set(ax,'YDir','normal');
    title(ax,label);
    xlabel(ax,unit);
    ylabel(ax,unit);
end
