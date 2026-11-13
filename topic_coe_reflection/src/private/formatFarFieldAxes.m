function formatFarFieldAxes(ax, cfg, style)
%formatFarFieldAxes 统一角度坐标，并标出雾度评价圆锥边界。
    axis(ax, 'image');
    ax.YDir = 'normal';
    limit = cfg.angleHalfRangeDeg;
    xlim(ax, [-limit, limit]); ylim(ax, [-limit, limit]);
    ax.XTick = linspace(-limit, limit, 5);
    ax.YTick = ax.XTick;
    xlabel(ax, 'θ_x / °'); ylabel(ax, 'θ_y / °');
    if ~style.showAnnotations
        return;
    end
    hold(ax, 'on');
    phi = linspace(0, 2*pi, 361);
    if style.showMetricCone
        plot(ax, cfg.metricAngleDeg*cos(phi), cfg.metricAngleDeg*sin(phi), ...
            '--', 'Color', [0.6, 0.85, 0.9], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
    obs = cfg.observation;
    if style.showObservationRegion
        plot(ax, obs.centerDeg(1) + obs.halfAngleDeg*cos(phi), obs.centerDeg(2) + obs.halfAngleDeg*sin(phi), ...
            ':', 'Color', [1, 0.7, 0.2], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
    if style.showObservationRegion && obs.innerHalfAngleDeg > 0
        plot(ax, obs.centerDeg(1) + obs.innerHalfAngleDeg*cos(phi), ...
            obs.centerDeg(2) + obs.innerHalfAngleDeg*sin(phi), ...
            ':', 'Color', [1, 0.7, 0.2], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
end
