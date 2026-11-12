function result = runCOEReflection(cfg)
%runCOEReflection COE black-state reflection, double-pass relief phase.
% cfg = coeReflectionDefaults(); edit cfg, then call here.
    arguments
        cfg (1,1) struct = coeReflectionDefaults()
    end
    assert(cfg.mode == "reflection", 'DisplayDiffraction:Mode', 'Use reflection configuration.');
    result = simulateDisplayDiffraction(cfg);
    baselineCfg = cfg;
    baselineCfg.heightMap = 0;
    baseline = simulateDisplayDiffraction(baselineCfg);
    result.baselineKernel = baseline.kernel;
    result.baselineMetrics = baseline.metrics;
    retention = NaN(height(result.metrics), 1);
    active = baseline.metrics.ZeroOrderFraction > 0;
    retention(active) = result.metrics.ZeroOrderFraction(active)./baseline.metrics.ZeroOrderFraction(active);
    result.metrics.ZeroOrderRetention = retention;
end
