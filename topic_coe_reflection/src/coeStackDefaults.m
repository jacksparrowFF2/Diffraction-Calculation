function stack = coeStackDefaults()
%coeStackDefaults 简化 Diamond COE 多层模型参数；所有材料参数均为可替换的示例值。
    arguments
    end
    stack.incidentN = 1;
    stack.anodeReflectance = 1;
    stack.anodePhaseDeg = 180; % 理想电反射边界，裸底 r=-1
    stack.bmExpansion = 5e-6;  % 阳极边界的 Euclidean 外扩距离，不是总宽度增加量
    stack.geometry.source = "diamond"; % diamond / mat / images
    stack.geometry.matFile = "coe_layers.mat";
    stack.geometry.anodeFile = "anode.png";
    stack.geometry.bmFile = ""; % 空表示按 bmExpansion 自动生成 BM 开口
    stack.geometry.bmImageIsOpening = true; % false：BM 图片白色表示实心 BM
    stack.geometry.cfFiles = ["CF_R.png", "CF_G.png", "CF_B.png"];
    stack.geometry.threshold = 0.5;
    stack.geometry.requireCFCoverage = true;
    stack.demo.anodeHalfDiagonal = [12, 10, 12]*1e-6; % R/G/B 菱形半对角线
    stack.demo.cfExpansion = 7e-6; % 大于 BM 外扩；过大会使不同 CF 交叠并报错
    stack.demo.colorTile = [1, 2; 2, 3]; % 2×2 演示晶胞 R/G/G/B，不代表厂家真实版图
    stack.materials.BM = material(1.8, 0.35, 1.0e-6);
    stack.materials.R = material(1.65, 0.03, 1.1e-6);
    stack.materials.G = material(1.60, 0.04, 1.0e-6);
    stack.materials.B = material(1.55, 0.05, 0.9e-6);
end

function value = material(n, k, thickness)
    value = struct('n', n, 'k', k, 'thickness', thickness, 'nkFile', "");
end
