function result = runCOEReflection(cfg)
%runCOEReflection COE black-state reflection, double-pass relief phase.
% cfg = coeReflectionDefaults(); edit cfg, then call here.
% 交互计算与结果保存请运行 COEReflection_main.m；此函数保留供批量计算调用。
% baseline 表示同一结构“不叠加额外高度相位”的状态，仍保留复反射系数自身的相位。
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
