function stack = coeLayerStackDefaults()
%coeLayerStackDefaults 从下到上的层级；相同 level 为同层互斥材料。
    stack.incidentN = 1;
    stack.exitN = 1; % 无基底有限膜：空气/膜堆/空气
    stack.backgroundMode = "absorber"; % COE 阳极外匹配吸收底；悬空膜用 medium
    stack.layerNames = ["Anode", "IJP", "BM", "CF_R", "CF_G", "CF_B"];
    stack.layerLevels = [1, 2, 3, 4, 4, 4];
    stack.geometry = struct('source', "images", 'matFile', "layers.mat", 'threshold', 0.5);
    for name = stack.layerNames
        stack.layers.(name) = coeLayerDefaults(name);
    end
    stack.layers.Anode.file = "anode";
    stack.layers.Anode.method = "reflectance";
    stack.layers.IJP.file = "IJP"; % 与其它图层一样读取同尺寸全白图片
    stack.layers.IJP.n = 1.5; stack.layers.IJP.k = 0;
    stack.layers.IJP.thickness = 10e-6;
    stack.layers.BM.n = 1.8; stack.layers.BM.k = 0.35; stack.layers.BM.thickness = 1e-6;
    stack.layers.CF_R.n = 1.65; stack.layers.CF_R.k = 0.03; stack.layers.CF_R.thickness = 1.1e-6;
    stack.layers.CF_G.n = 1.60; stack.layers.CF_G.k = 0.04; stack.layers.CF_G.thickness = 1e-6;
    stack.layers.CF_B.n = 1.55; stack.layers.CF_B.k = 0.05; stack.layers.CF_B.thickness = 0.9e-6;
end
