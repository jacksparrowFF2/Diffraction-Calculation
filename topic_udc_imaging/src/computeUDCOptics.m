function result = computeUDCOptics(cfg)
%computeUDCOptics Calculate transmission kernels and MTF without a scene image.
    arguments
        cfg (1,1) struct = udcImagingDefaults()
    end
    assert(cfg.mode == "udc", 'DisplayDiffraction:Mode', 'Use udc configuration.');
    result = simulateDisplayDiffraction(cfg);
end
