function [result, baseline, heightPhaseEnabled] = coreSimulateHeightStates(cfg, progressFcn)
%coreSimulateHeightStates 额外高度全零时只传播一次，兼容保留 baseline 数组字段。
    if nargin < 2, progressFcn = @(fraction,message) []; end
    heightPhaseEnabled = any(cfg.heightMap(:) ~= 0);
    stateCount = 1+double(heightPhaseEnabled);
    stateLabel = "反射结果（额外高度 0 nm）";
    if heightPhaseEnabled, stateLabel = "叠加高度相位"; end
    result = coreSimulateDiffraction(cfg, ...
        @(fraction,message) progressFcn(fraction/stateCount, stateLabel+"｜"+message));
    if heightPhaseEnabled
        baselineCfg = cfg;
        baselineCfg.heightMap = 0;
        baseline = coreSimulateDiffraction(baselineCfg, ...
            @(fraction,message) progressFcn((1+fraction)/stateCount, "不叠加高度相位｜"+message));
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
