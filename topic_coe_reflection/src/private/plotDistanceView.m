function [fig,displayData] = plotDistanceView(view,spectrum,colorOptions,rgbOptions,isColorimetric)
%plotDistanceView 独立显示几何观察角图，避免与远场接收面坐标混用。
    baseline=view;
    baseline.spectralEnergy=view.baselineSpectralEnergy;
    if isColorimetric
        [~,color,mapping]=coreSynthesizeColor(baseline,view,spectrum,colorOptions,view.heightPhaseEnabled);
        rgb=color.displaySRGB;
        label="光谱积分sRGB";
    else
        assert(size(view.channelSignal,3)==3,'COEReflection:DistanceViewChannels', ...
            '距离视图预览要求RGB三通道或光谱色度积分。');
        [~,rgb,mapping]=plotComposeReflectionRgb(view.baselineChannelSignal/view.binSolidAngle, ...
            view.channelSignal/view.binSolidAngle,rgbOptions,view.heightPhaseEnabled);
        label="三波长RGB代理";
    end
    fig=figure('Name','COE 距离缩放观察视图','Color','w');
    imagesc(rad2deg(view.coordinate),rad2deg(view.coordinate),rgb);
    axis image xy; xlabel('观察角 theta_x / °'); ylabel('观察角 theta_y / °');
    title({sprintf('距离缩放视图 · %d×%d Pitch · a=%.1f cm · b=%.1f cm · s=%.4f', ...
        view.config.periods(1),view.config.periods(2),100*view.sourceDistance,100*view.observerDistance,view.angularScale), ...
        label+sprintf('；源%.1f mm，入瞳%.1f mm；%s，非完整眼成像', ...
        view.sourceDiameter*1e3,view.pupilDiameter*1e3,view.kernels.pupilModel)});
    displayData=struct('rgb',rgb,'mapping',mapping);
end
