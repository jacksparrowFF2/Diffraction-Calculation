function styleCoeFigure(fig, style)
%styleCoeFigure 统一中文字体、坐标轴和色标样式，仅作用于本次生成的图。
    set(findall(fig, '-property', 'FontName'), 'FontName', style.fontName);
    set(findall(fig, '-property', 'FontSize'), 'FontSize', style.fontSize);
    set(findall(fig, 'Type', 'axes'), 'LineWidth', style.lineWidth, ...
        'TickDir', 'out', 'Box', 'on', 'TickLabelInterpreter', 'none');
    set(findall(fig, 'Type', 'colorbar'), 'LineWidth', style.lineWidth, 'TickDirection', 'out');
    set(findall(fig, 'Type', 'legend'), 'Box', 'off');
    set(findall(fig, 'Tag', 'COESpotLabel'), 'FontSize', 8.5);
end
