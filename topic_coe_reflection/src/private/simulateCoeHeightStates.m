function [result, baseline, heightPhaseEnabled] = simulateCoeHeightStates(cfg)
%simulateCoeHeightStates 额外高度全零时只传播一次，兼容保留 baseline 数组字段。
    result = simulateDisplayDiffraction(cfg);
    heightPhaseEnabled = any(cfg.heightMap(:) ~= 0);
    if heightPhaseEnabled
        baselineCfg = cfg;
        baselineCfg.heightMap = 0;
        baseline = simulateDisplayDiffraction(baselineCfg);
    else
        baseline = result; % MATLAB 写时复制；复用结果，不再执行传播或膜厚计算
    end
    result.heightPhaseEnabled = heightPhaseEnabled;
    result.propagationRuns = 1 + double(heightPhaseEnabled);
    result.baselineKernel = baseline.kernel;
    result.baselineMetrics = baseline.metrics;
    retention = NaN(height(result.metrics),1);
    active = baseline.metrics.ZeroOrderFraction > 0;
    if heightPhaseEnabled
        retention(active) = result.metrics.ZeroOrderFraction(active)./baseline.metrics.ZeroOrderFraction(active);
    else
        retention(active) = 1;
    end
    result.metrics.ZeroOrderRetention = retention;
end
