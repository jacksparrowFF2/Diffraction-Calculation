function reportCoeSampling(cfg, model)
%reportCoeSampling 说明图像网格与厚膜光谱采样；不为消除提示而修改物理参数。
    extent=max(cfg.periods)*cfg.pitch;
    step=2*cfg.angleHalfRangeDeg/(cfg.outputSize-1);
    width=rad2deg(min(cfg.wavelengths)/extent);
    fprintf('角域 ±%.3f°，%d×%d 网格，间距 %.5f°；最短波长 λ/D≈%.5f°，约 %.2f 个显示间隔。\n', ...
        cfg.angleHalfRangeDeg,cfg.outputSize,cfg.outputSize,step,width,width/step);
    fprintf('输出像素内部仍自动细分积分；局部放大只裁切显示，不补充真实角度分辨率。\n');
    if numel(cfg.wavelengths)>1
        materialsTable=model.materialTable;
        opticalDepth=0;
        for level=unique(materialsTable.Level).'
            selected=materialsTable.Level==level & materialsTable.Method=="nk";
            if any(selected), opticalDepth=opticalDepth+max(materialsTable.n(selected).*materialsTable.ThicknessM(selected)); end
        end
        phaseStep=4*pi*opticalDepth*max(abs(diff(1./cfg.wavelengths)));
        if phaseStep>pi/2
            suggested=min(cfg.wavelengths)^2/(8*opticalDepth)*1e9;
            warning('COEReflection:ThickFilmSpectrum', ...
                '厚膜的光谱干涉可能欠采样（最大往返相位步进约 %.2f rad）。当前可作预览；材料/色度定量比较建议波长步长约 %.3g nm 或更小，并检查收敛。', ...
                phaseStep,suggested);
        end
    end
end
