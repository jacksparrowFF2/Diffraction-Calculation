function source = coeSourceDefaults()
%coeSourceDefaults 光源空间模型；光谱独立设置，入射辐照度为物面平均值。
    source.type = "plane"; % plane / point / laser
    source.centerXY = [0,0]; % 点源中心或高斯光轴在物面的投影，m
    source.distance = 0.1; % 点源到物面法向距离，m
    source.diameter = 0; % 点源发光圆盘直径，0=理想点；m
    source.samplesPerAxis = 3; % 有限圆盘确定性采样；各点互不相干
    source.beamHalfAngleDeg = 30; % 点源辐射强度降至轴向一半的半角
    source.laserSpecification = "waist"; % waist 或 divergence，避免独立指定矛盾参数
    source.waistRadius = 100e-6; % 基模高斯束腰1/e²强度半径，m
    source.waistDistance = 0; % 沿传播方向从束腰到样品的有符号距离，m
    source.divergenceHalfAngleDeg = 0.1; % 激光1/e²强度远场发散半角
end
