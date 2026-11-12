function heightMap = createCOEHeightExample(cfg, example)
%createCOEHeightExample 生成可复现的二值高度条纹，单位 m，不改变材料反射率分布。
% 条纹从计算矩阵首行/列开始排列；周期和高区宽度必须能被当前采样精确表示。
    arguments
        cfg (1,1) struct
        example (1,1) struct
    end
    validateattributes(example.height, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
    validateattributes(example.period, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
    validateattributes(example.dutyCycle, {'numeric'}, {'scalar', 'real', '>', 0, '<', 1});
    mustBeMember(example.direction, ["x", "y"]);
    dx = cfg.pitch/cfg.samplesPerPeriod;
    periodSamples = example.period/dx;
    highSamples = periodSamples*example.dutyCycle;
    assert(periodSamples >= 4 && abs(periodSamples - round(periodSamples)) < 1e-9 ...
        && abs(highSamples - round(highSamples)) < 1e-9, ...
        'COEReflection:HeightSampling', ...
        '高度条纹周期至少需要 4 个采样点，周期及高区宽度均须为采样间距的整数倍。');
    periodSamples = round(periodSamples);
    highSamples = round(highSamples);
    gridSize = [cfg.periods(2), cfg.periods(1)]*cfg.samplesPerPeriod;
    if example.direction == "x"
        highRegion = mod(0:gridSize(2)-1, periodSamples) < highSamples;
        heightMap = repmat(double(highRegion)*example.height, gridSize(1), 1);
    else
        highRegion = mod((0:gridSize(1)-1).', periodSamples) < highSamples;
        heightMap = repmat(double(highRegion)*example.height, 1, gridSize(2));
    end
end
