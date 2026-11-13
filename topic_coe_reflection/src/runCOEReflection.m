function result = runCOEReflection(cfg, stack, inputDir)
%runCOEReflection COE black-state reflection, double-pass relief phase.
% runCOEReflection(cfg, stack, inputDir) 按统一图层列表读取图片并计算。
% 仅传 cfg 时保留历史调用约定，不自动读取项目 config.txt。
% 交互计算与结果保存请运行 COEReflection_main.m；此函数保留供批量计算调用。
% baseline 表示同一结构“不叠加额外高度相位”的状态，仍保留复反射系数自身的相位。
    arguments
        cfg (1,1) struct = coeReflectionDefaults()
        stack (1,1) struct = struct()
        inputDir (1,1) string = ""
    end
    assert(cfg.mode == "reflection", 'DisplayDiffraction:Mode', 'Use reflection configuration.');
    if ~isempty(fieldnames(stack))
        [cfg.coeStack, geometry] = prepareCoeLayerStack(cfg, stack, inputDir);
        cfg.mask = double(geometry.unionMask);
    end
    result = simulateCoeHeightStates(cfg);
    if ~isempty(fieldnames(stack))
        result.layerGeometry = geometry;
    end
end
