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
    options = optics.quality.options;
    floorDb = [];
    if isfield(options,'rgbDisplayFloorDb'), floorDb = options.rgbDisplayFloorDb; end
    ax = nexttile(layout);
    reference = max(max(optics.clearKernel,[],'all'),realmin);
    if isempty(floorDb)
        preview = min(log1p(1000*optics.kernel/reference)/log1p(1000),1);
        rgbTitle = 'RGB 衍射强度（对数显示）';
    else
        preview = rgbInDb(optics.kernel,reference,floorDb);
        rgbTitle = sprintf('RGB 衍射强度（%g 至 0 dB）',floorDb);
    end
    drawMap(ax,coordinate,coordinate,preview,rgbTitle,'传感器 (μm)');

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
    selected = abs(coordinate)<=options.hdrDisplayHalfWidthUm+8*eps(max(abs(coordinate)));
    displayCoordinate = coordinate(selected);
    clearPoint = point.clear(selected,selected,:);
    maskedPoint = point.withMask(selected,selected,:);
    heading = sprintf('HDR 点光源 · 输入 %.3g · 共同曝光',options.hdrPointSignal);
    if isempty(floorDb)
        reference = max([1,max(point.clear,[],'all'),max(point.withMask,[],'all')]);
        displayScale = log1p(20*reference);
        clearPreview = log1p(20*clearPoint)/displayScale;
        maskedPreview = log1p(20*maskedPoint)/displayScale;
        title(layout,string(heading)+"与对数显示");
    else
        reference = max(max(point.clear,[],'all'),realmin);
        clearPreview = rgbInDb(clearPoint,reference,floorDb);
        maskedPreview = rgbInDb(maskedPoint,reference,floorDb);
        title(layout,{heading,sprintf('%g 至 0 dB · 0 dB 为无遮挡峰值',floorDb)});
    end
    drawMap(nexttile(layout),displayCoordinate,displayCoordinate,clearPreview,'无遮挡','传感器 (μm)');
    drawMap(nexttile(layout),displayCoordinate,displayCoordinate,maskedPreview,'含 UDC 掩膜','传感器 (μm)');
    exportgraphics(layout,fullfile(outputDir,'udc_hdr_point.png'),'Resolution',160);
end

function rgb = rgbInDb(energy,reference,floorDb)
    % 按每个像元最强通道映射亮度，保留该像元原来的 RGB 比例。
    peak = max(max(energy,[],3),realmin);
    db = 10*(log10(peak)-log10(reference));
    level = min(max(1-db/floorDb,0),1);
    rgb = (energy./peak).*level;
end

function drawMap(ax,x,y,data,label,unit)
    imagesc(ax,x,y,data);
    axis(ax,'image');
    set(ax,'YDir','normal');
    title(ax,label);
    xlabel(ax,unit);
    ylabel(ax,unit);
end
