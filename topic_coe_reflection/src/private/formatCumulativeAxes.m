function formatCumulativeAxes(observation, labels, panelTitle, style)
%formatCumulativeAxes 半径内累计值始终以入射量为分母，不强制将曲线末端归一化到 100%。
    if style.showAnnotations && style.showObservationRegion
        xline(observation.halfAngleDeg, ':', '观察外半角', 'HandleVisibility', 'off');
    end
    if style.showAnnotations && style.showObservationRegion && observation.innerHalfAngleDeg > 0
        xline(observation.innerHalfAngleDeg, ':', '观察内半角', 'HandleVisibility', 'off');
    end
    xlabel('相对观察中心的累计半角 / °'); ylabel('累计绝对反射率 / %');
    title(panelTitle); legend(labels, 'Location', 'best'); grid on;
    xlim([0, max(xlim)]);
end
