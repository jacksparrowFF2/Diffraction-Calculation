function angular = computeCOEAngularMetrics(result, observation, spectrum)
%computeCOEAngularMetrics 对观察圆锥/环带积分，计算以入射能量为分母的累计反射率。
% 角域按输出 bin 中心判定；角边界附近应检查 outputSize 和 FFT 采样收敛。
% energyWeights 为辐射度权重，photopicWeights 为明视觉权重，二者均归一化为单位入射量。
    arguments
        result (1,1) struct
        observation (1,1) struct
        spectrum (1,1) struct
    end
    validateattributes(observation.centerDeg, {'numeric'}, {'real', 'finite', 'size', [1, 2]});
    validateattributes(observation.halfAngleDeg, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
    validateattributes(observation.innerHalfAngleDeg, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'nonnegative', '<', observation.halfAngleDeg});
    validateattributes(observation.curveSamples, {'numeric'}, {'scalar', 'integer', '>=', 2});
    angleDeg = rad2deg(result.coordinate);
    assert(all(abs(angleDeg + fliplr(angleDeg)) < 1e-9), ...
        'COEReflection:AngleGrid', '累计角能量要求对称的反射角网格。');
    limit = max(angleDeg);
    maxRadius = min(limit - abs(observation.centerDeg));
    assert(observation.halfAngleDeg <= maxRadius, 'COEReflection:ObservationExtent', ...
        '观察角域超出输出窗口；请移动中心、减小半角或扩大观察窗口。');
    [xx, yy] = meshgrid(angleDeg);
    radius = hypot(xx - observation.centerDeg(1), yy - observation.centerDeg(2));
    mask = radius <= observation.halfAngleDeg;
    if observation.innerHalfAngleDeg > 0
        mask = mask & radius > observation.innerHalfAngleDeg;
    end
    signals = reshape(result.spectralEnergy, [], numel(result.config.wavelengths));
    nLambda = size(signals, 2);
    energyWeights = spectrum.energyWeights(:);
    photopicWeights = spectrum.photopicWeights(:);
    validateattributes(energyWeights, {'numeric'}, {'numel', nLambda, 'finite', 'nonnegative'});
    validateattributes(photopicWeights, {'numeric'}, {'numel', nLambda, 'finite', 'nonnegative'});
    assert(abs(sum(energyWeights) - 1) < 1e-10 && abs(sum(photopicWeights) - 1) < 1e-10, ...
        'COEReflection:AngularWeights', '积分权重必须归一化为单位入射量。');
    perLambda = sum(signals(mask(:), :), 1).';
    angular.observation = observation;
    angular.mask = mask;
    angular.spectralViewReflectance = perLambda;
    angular.radiantViewReflectance = perLambda.'*energyWeights;
    angular.photopicViewReflectance = perLambda.'*photopicWeights;
    angular.radiantTotalReflectance = result.metrics.Throughput.'*energyWeights;
    angular.photopicTotalReflectance = result.metrics.Throughput.'*photopicWeights;
    angular.radiantCapturedReflectance = result.metrics.CapturedThroughput.'*energyWeights;
    angular.photopicCapturedReflectance = result.metrics.CapturedThroughput.'*photopicWeights;
    % 累计曲线从观察中心积分到半径 alpha，不减去观察环带的内孔。
    radii = linspace(0, maxRadius, observation.curveSamples);
    % 精确保留用户设定边界，移除浮点舍入造成的近重复半径。
    nearBoundary = abs(radii - observation.innerHalfAngleDeg) < 1e-10 ...
        | abs(radii - observation.halfAngleDeg) < 1e-10;
    radii(nearBoundary) = [];
    radii = unique([radii, observation.innerHalfAngleDeg, observation.halfAngleDeg]);
    [sortedRadius, order] = sort(radius(:));
    counts = zeros(size(radii));
    for j = 1:numel(radii)
        if radii(j) > 0
            counts(j) = nnz(sortedRadius <= radii(j));
        end
    end
    cumulative = zeros(numel(radii), nLambda);
    active = counts > 0;
    for k = 1:nLambda
        running = cumsum(signals(order, k));
        cumulative(active, k) = running(counts(active));
    end
    angular.radiusDeg = radii(:);
    angular.spectralCumulativeReflectance = cumulative;
    angular.radiantCumulativeReflectance = cumulative*energyWeights;
    angular.photopicCumulativeReflectance = cumulative*photopicWeights;
end
