function markDiffractionSpots(ax, spots, groupName, wavelength, style)
%markDiffractionSpots 标出共同积分圆、编号与零级排除圆；不使用显示增强结果找峰。
    if ~style.showAnnotations
        return;
    end
    groupIndex = find(string({spots.groups.name}) == groupName);
    if groupName == "wavelength"
        groupIndex = groupIndex(abs([spots.groups(groupIndex).wavelengthNm] - wavelength) < 1e-6);
    end
    if isempty(groupIndex)
        return;
    end
    centers = spots.groups(groupIndex(1)).centersDeg;
    phi = linspace(0, 2*pi, 181);
    hold(ax, 'on');
    if style.showZeroOrderExclusion
        plot(ax, spots.options.zeroOrderExclusionDeg*cos(phi), ...
            spots.options.zeroOrderExclusionDeg*sin(phi), ':', 'Color', [1, 0.35, 0.35], ...
            'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
    for j = 1:size(centers, 1)
        if style.showSpotRegions
            plot(ax, centers(j, 1) + spots.options.integrationHalfAngleDeg*cos(phi), ...
                centers(j, 2) + spots.options.integrationHalfAngleDeg*sin(phi), ...
                '-', 'Color', [1, 0.35, 0.8], 'LineWidth', 0.8, 'HandleVisibility', 'off');
        end
        if ~style.showSpotLabels
            continue;
        end
        limitsX = xlim(ax); limitsY = ylim(ax);
        if centers(j, 1) >= limitsX(1) && centers(j, 1) <= limitsX(2) ...
                && centers(j, 2) >= limitsY(1) && centers(j, 2) <= limitsY(2)
            direction = centers(j, :)/max(norm(centers(j, :)), eps);
            span = max(diff(limitsX), diff(limitsY));
            offset = max(1.3*spots.options.integrationHalfAngleDeg, 0.065*span);
            labelPosition = centers(j, :) + offset*direction;
            labelPosition(1) = min(max(labelPosition(1), limitsX(1) + 0.04*span), limitsX(2) - 0.04*span);
            labelPosition(2) = min(max(labelPosition(2), limitsY(1) + 0.04*span), limitsY(2) - 0.04*span);
            if offset > 1.5*spots.options.integrationHalfAngleDeg
                start = centers(j, :) + spots.options.integrationHalfAngleDeg*direction;
                plot(ax, [start(1), labelPosition(1)], [start(2), labelPosition(2)], ...
                    '-', 'Color', [0.8, 0.3, 0.65], 'LineWidth', 0.5, 'HandleVisibility', 'off');
            end
            text(ax, labelPosition(1), labelPosition(2), "S" + j, ...
                'Color', [1, 0.75, 0.9], 'HorizontalAlignment', 'center', ...
                'FontSize', 8.5, 'Tag', 'COESpotLabel', 'Clipping', 'on');
        end
    end
end
