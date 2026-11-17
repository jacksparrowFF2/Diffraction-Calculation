function result = runCOEReflection(cfg, stack, inputDir)
%runCOEReflection COE black-state reflection, double-pass relief phase.
% runCOEReflection(cfg, stack, inputDir) 按统一图层列表读取图片并计算。
% 仅传 cfg 时保留历史调用约定，不自动读取项目 config.txt。
% 交互计算与结果保存请运行 COEReflection_main.m；此函数保留供批量计算调用。
% baseline 表示同一结构“不叠加额外高度相位”的状态，仍保留复反射系数自身的相位。
% cfg.observer.enabled=true时，额外返回独立几何观察视图result.distanceView。
% source.distance为光源距离，observer.distance为眼睛距离；不替代propagationDistance。
    arguments
        cfg (1,1) struct = coeReflectionDefaults()
        stack (1,1) struct = struct()
        inputDir (1,1) string = ""
    end
    assert(cfg.mode == "reflection", 'DisplayDiffraction:Mode', 'Use reflection configuration.');
    if ~isempty(fieldnames(stack))
        [cfg.coeStack, geometry] = corePrepareLayerStack(cfg, stack, inputDir);
        cfg.mask = double(geometry.unionMask);
    end
    if isfield(cfg,'observer') && cfg.observer.enabled && isfield(cfg.observer,'model') && cfg.observer.model=="coherentEye"
        result=struct('eyeImage',coreEyeImage(cfg));
        if ~isempty(fieldnames(stack)), result.layerGeometry=geometry; end
        return
    end
    result = coreSimulateHeightStates(cfg);
    if isfield(cfg,'observer') && cfg.observer.enabled
        result.distanceView = coreDistanceView(cfg);
    end
    if ~isempty(fieldnames(stack))
        result.layerGeometry = geometry;
    end
end
