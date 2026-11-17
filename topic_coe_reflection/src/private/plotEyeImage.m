function [fig,displayData]=plotEyeImage(view,rgbOptions)
%plotEyeImage 像面坐标按retinaDistance换算观察角；全视野与裁看共用曝光。
    angle=rad2deg(view.coordinate); cfg=view.config;
    referenceHalf=0.25;
    if isfield(cfg.observer,'displayReferenceHalfAngleDeg'), referenceHalf=cfg.observer.displayReferenceHalfAngleDeg; end
    validateattributes(referenceHalf,{'numeric'},{'scalar','real','finite','nonnegative'});
    if referenceHalf>0
        selected=abs(angle)<=referenceHalf+1e-12;
        assert(any(selected),'COEReflection:EyeDisplayReference','显示参考ROI无采样点。');
        [~,~,referenceMapping]=plotComposeReflectionRgb( ...
            view.baselineChannelSignal(selected,selected,:)/view.binSolidAngle, ...
            view.channelSignal(selected,selected,:)/view.binSolidAngle,rgbOptions,view.heightPhaseEnabled);
        rgbOptions.fixedReferenceSignal=referenceMapping.referenceSignal;
    end
    [~,rgb,mapping]=plotComposeReflectionRgb(view.baselineChannelSignal/view.binSolidAngle, ...
        view.channelSignal/view.binSolidAngle,rgbOptions,view.heightPhaseEnabled);
    mapping.referenceHalfAngleDeg=referenceHalf;
    fig=figure('Name','逐源点复场人眼成像','Color','w','Position',[40,70,1450,740]);
    tiles=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
    title(tiles,sprintf('逐源点复场→圆瞳/薄透镜→像面强度和 · %d×%d Pitch',cfg.periods));
    focusLabel=cfg.observer.focusMode;
    if focusLabel=="virtualSource", focusLabel="镜中光源"; elseif focusLabel=="screen", focusLabel="屏幕"; end
    subtitle(tiles,sprintf('源%.1f mm · 入瞳%.1f mm · a/b=%.0f/%.0f cm · 对焦%s，%.0f cm · %s RGB代理\n共同曝光参考±%.2f°；无空间亮度遮罩', ...
        cfg.source.diameter*1e3,cfg.observer.pupilDiameter*1e3,cfg.source.distance*100, ...
        cfg.observer.distance*100,focusLabel,view.focusDistance*100,rgbOptions.mode,referenceHalf));
    for j=1:2
        nexttile; imagesc(angle,angle,rgb); axis image xy;
        xlabel('像面位置 / 像距 = 观察角 theta_x / °'); ylabel('观察角 theta_y / °');
        if j==1
            title(sprintf('完整视野 ±%.2f°',max(angle)));
        else
            half=cfg.observer.imageZoomHalfAngleDeg;
            xlim([-half,half]); ylim([-half,half]); title(sprintf('±%.2f°裁看；相同曝光和原始分辨率',half));
        end
    end
    displayData=struct('rgb',rgb,'mapping',mapping);
end
