function [comparisonFigure,receiverFigure,view]=plotReceiverDistanceViews(result,rgb,distances)
%plotReceiverDistanceViews 同一远场角分布映射到X=z*theta；毫米尺度共轴比较。
% RGB归一形状显示；物理辐照度按(z_ref/z)^2缩放并单独记录，窗口能量不重标定。
    validateattributes(distances,{'numeric'},{'row','finite','positive','nonempty','increasing'});
    assert(result.config.mode=="reflection",'COEReflection:ReceiverViewMode','接收尺寸视图要求反射远场角坐标。');
    assert(isequal(size(rgb,1:2),size(result.spectralEnergy,1:2)) && size(rgb,3)==3, ...
        'COEReflection:ReceiverViewRGB','接收尺寸视图要求同网格RGB图。');
    theta=result.coordinate; dt=theta(2)-theta(1);
    referenceDistance=result.config.propagationDistance;
    coordinateMm=distances(:)*theta(:).'*1e3;
    irradianceScale=(referenceDistance./distances).^2;
    [~,green]=min(abs(result.config.wavelengths-550e-9));
    centre=(numel(theta)+1)/2;
    centreIrradiance=result.spectralIrradiance(centre,centre,green)*irradianceScale;
    greenPower=result.incidentPower*result.metrics.CapturedThroughput(green);
    nominalOrder=distances*result.config.wavelengths(green)/result.config.pitch*1e3;
    halfWidth=max(abs(coordinateMm),[],2);
    summary=table(distances(:),halfWidth,2*halfWidth,nominalOrder(:), ...
        centreIrradiance(:),repmat(greenPower,numel(distances),1), ...
        'VariableNames',{'ReceiverDistanceM','HalfWidthMm','FullWidthMm', ...
        'NominalGreenFirstOrderMm','GreenCentreIrradianceWm2','GreenWindowPowerW'});
    view=struct('distances',distances,'coordinateMm',coordinateMm, ...
        'referenceDistanceM',referenceDistance,'irradianceScaleFromReference',irradianceScale, ...
        'pixelAreaM2',(distances*dt).^2,'summary',summary, ...
        'coordinateConvention',"X=z*theta, paraxial far-field projection", ...
        'displayConvention',"same normalized RGB shape; absolute irradiance changes are in summary", ...
        'assumptions',"固定平面波照明；同一远场角分布的坐标换算，未验证这些z处的近场传播");
    count=numel(distances); limits=[-max(halfWidth),max(halfWidth)];
    comparisonFigure=figure('Name','接收距离与实际图样尺寸','Color','w','Position',[50,80,1450,570]);
    tiles=tiledlayout(1,count,'TileSpacing','compact','Padding','compact');
    title(tiles,sprintf('平面波 · %d×%d Pitch · X=z theta · 所有面板共用毫米坐标',result.config.periods));
    subtitle(tiles,'同一远场角图的尺寸映射；RGB归一显示以比较形状，绝对辐照度见CSV');
    for j=1:count
        ax=nexttile; imagesc(coordinateMm(j,:),coordinateMm(j,:),rgb); axis image xy;
        xlim(limits); ylim(limits); set(ax,'Color','k');
        xlabel('接收面 X / mm'); ylabel('接收面 Y / mm');
        title(sprintf('z=%.0f cm · 角域窗口全宽%.2f mm',distances(j)*100,2*halfWidth(j)));
    end
    receiverFigure=figure('Name','接收面实际尺寸RGB','Color','w','Position',[70,60,950,940]);
    imagesc(coordinateMm(1,:),coordinateMm(1,:),rgb); axis image xy;
    xlabel('接收面 X / mm'); ylabel('接收面 Y / mm');
    title(sprintf('平面波 · %d×%d Pitch · 光栅到接收面 z=%.0f cm',result.config.periods,distances(1)*100));
    subtitle(sprintf('X=z theta · 角域±%.1f°对应±%.2f mm · RGB形状显示',rad2deg(max(theta)),halfWidth(1)));
end
