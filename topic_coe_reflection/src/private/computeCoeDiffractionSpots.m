function spots = computeCoeDiffractionSpots(withoutResult, withResult, spectrum, options, compareStates)
%computeCoeDiffractionSpots 在原始能量图上检测非零级亮斑，以固定圆域积分绝对反射率。
% 两状态的逐点最大能量图仅用于共同定位和排序，实际能量始终分别从各自原始数组积分。
% 相同组内的积分圆不重叠、不进入零级排除圆，也不越过输出窗口，不能把不同组的数值相加。
    arguments
        withoutResult (1,1) struct
        withResult (1,1) struct
        spectrum (1,1) struct
        options (1,1) struct
        compareStates (1,1) logical = true
    end
    validateattributes(options.zeroOrderExclusionDeg, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'});
    validateattributes(options.integrationHalfAngleDeg, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
    validateattributes(options.minSeparationDeg, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'});
    validateattributes(options.maxSpots, {'numeric'}, {'scalar', 'integer', 'positive'});
    validateattributes(options.minRelativePeak, {'numeric'}, {'scalar', 'real', '>=', 0, '<=', 1});
    if ~isempty(options.wavelengthNm)
        validateattributes(options.wavelengthNm, {'numeric'}, {'vector', 'real', 'finite', 'positive'});
    end
    assert(isequal(withoutResult.coordinate, withResult.coordinate) ...
        && isequal(withoutResult.config.wavelengths, withResult.config.wavelengths), ...
        'COEReflection:SpotGrid', '两状态必须使用相同角度与波长网格。');
    coordinate = rad2deg(withResult.coordinate);
    step = coordinate(2) - coordinate(1);
    assert(options.integrationHalfAngleDeg >= step, 'COEReflection:SpotSampling', ...
        '衍射斑积分半角不得小于一个输出网格间距（当前 %.4f°）。', step);
    assert(options.zeroOrderExclusionDeg + 2*options.integrationHalfAngleDeg < max(coordinate), ...
        'COEReflection:SpotExtent', '零级排除范围和积分半角过大，输出窗口中没有完整的候选区域。');
    nLambda = numel(withResult.config.wavelengths);
    a = reshape(withoutResult.spectralEnergy, [], nLambda);
    b = reshape(withResult.spectralEnergy, [], nLambda);
    if isempty(options.wavelengthNm)
        selected = 1:nLambda;
    else
        selected = zeros(size(options.wavelengthNm));
        for j = 1:numel(selected)
            [~, selected(j)] = min(abs(withResult.config.wavelengths*1e9 - options.wavelengthNm(j)));
        end
        selected = unique(selected, 'stable');
    end
    weights = [spectrum.energyWeights(:), spectrum.photopicWeights(:)];
    assert(size(weights, 1) == nLambda && all(weights >= 0, 'all') ...
        && all(abs(sum(weights, 1) - 1) < 1e-10), ...
        'COEReflection:SpotWeights', '光谱积分权重必须非负且归一化为单位入射量。');
    [xx, yy] = meshgrid(coordinate);
    groups = cell(1, numel(selected) + 2);
    for j = 1:numel(selected)
        k = selected(j);
        groups{j} = analyseGroup(reshape(a(:, k), size(xx)), reshape(b(:, k), size(xx)), ...
            xx, yy, options, "wavelength", withResult.config.wavelengths(k)*1e9, compareStates);
    end
    if ~compareStates, b = a; end
    names = ["radiant", "photopic"];
    for j = 1:2
        firstMap = reshape(a*weights(:, j), size(xx));
        secondMap = firstMap;
        if compareStates, secondMap = reshape(b*weights(:,j),size(xx)); end
        groups{numel(selected)+j} = analyseGroup(firstMap,secondMap,xx,yy,options,names(j),NaN,compareStates);
    end
    spots.options = options;
    spots.effectiveMinSeparationDeg = max(options.minSeparationDeg, 2*options.integrationHalfAngleDeg);
    spots.groups = [groups{:}];
    tables = arrayfun(@(group) group.metrics, spots.groups, 'UniformOutput', false);
    spots.metrics = vertcat(tables{:});
end

function group = analyseGroup(a, b, xx, yy, options, name, wavelength, compareStates)
    reference = max(a, b);
    n = size(reference, 1);
    padded = -Inf(n + 2);
    padded(2:end-1, 2:end-1) = reference;
    peak = true(size(reference));
    strict = false(size(reference));
    for dy = -1:1
        for dx = -1:1
            if dx ~= 0 || dy ~= 0
                neighbour = padded((2:n+1) + dy, (2:n+1) + dx);
                peak = peak & reference >= neighbour;
                strict = strict | reference > neighbour;
            end
        end
    end
    halfAngle = options.integrationHalfAngleDeg;
    eligible = hypot(xx, yy) >= options.zeroOrderExclusionDeg + halfAngle ...
        & abs(xx) + halfAngle <= max(xx(:)) & abs(yy) + halfAngle <= max(yy(:));
    strongest = max(reference(eligible), [], 'all');
    if isempty(strongest)
        strongest = 0;
    end
    candidates = find(peak & strict & eligible & reference > 0 ...
        & reference >= options.minRelativePeak*strongest);
    % 先按候选积分能量排序，同分时按线性索引固定顺序，避免同一对比中标签随机变化。
    score = zeros(numel(candidates), 1);
    for j = 1:numel(candidates)
        circle = hypot(xx - xx(candidates(j)), yy - yy(candidates(j))) <= halfAngle;
        score(j) = sum(reference(circle));
    end
    [~, ranking] = sortrows([-score, double(candidates)], [1, 2]);
    chosen = zeros(options.maxSpots, 1);
    count = 0;
    separation = max(options.minSeparationDeg, 2*halfAngle);
    for index = ranking(:).'
        candidate = candidates(index);
        if count == 0 || all(hypot(xx(chosen(1:count)) - xx(candidate), ...
                yy(chosen(1:count)) - yy(candidate)) > separation)
            count = count + 1;
            chosen(count) = candidate;
            if count == options.maxSpots
                break;
            end
        end
    end
    chosen = chosen(1:count);
    centers = [xx(chosen), yy(chosen)];
    stateCount = 1+double(compareStates);
    energyPercent = zeros(stateCount*count, 1);
    peakPosition = NaN(stateCount*count, 2);
    maps = {a, b};
    occupied = false(size(a));
    for j = 1:count
        % 按距离去重；精确落在圆周上的共享 bin 只归前一个 ROI，避免重复计数。
        circle = hypot(xx - centers(j, 1), yy - centers(j, 2)) <= halfAngle & ~occupied;
        occupied = occupied | circle;
        for state = 1:stateCount
            row = stateCount*(j-1) + state;
            values = maps{state};
            energyPercent(row) = 100*sum(values(circle));
            masked = values;
            masked(~circle) = -Inf;
            [maximum, location] = max(masked(:));
            if maximum > 0
                peakPosition(row, :) = [xx(location), yy(location)];
            end
        end
    end
    id = repelem((1:count).', stateCount);
    stateIndex = repmat((1:stateCount).', count, 1);
    stateLabels = ["不叠加高度相位"; "叠加高度相位"];
    if ~compareStates, stateLabels = "反射结果（额外高度 0 nm）"; end
    centerRows = repelem(centers, stateCount, 1);
    metrics = table(repmat(name, stateCount*count, 1), repmat(wavelength, stateCount*count, 1), id, stateIndex, ...
        stateLabels(stateIndex), centerRows(:, 1), centerRows(:, 2), hypot(centerRows(:, 1), centerRows(:, 2)), ...
        peakPosition(:, 1), peakPosition(:, 2), repmat(halfAngle, stateCount*count, 1), ...
        repmat(options.zeroOrderExclusionDeg, stateCount*count, 1), energyPercent, ...
        'VariableNames', {'Group', 'WavelengthNm', 'SpotID', 'StateIndex', 'HeightPhaseState', ...
        'CenterThetaXDeg', 'CenterThetaYDeg', 'CenterRadiusDeg', 'PeakThetaXDeg', 'PeakThetaYDeg', ...
        'IntegrationHalfAngleDeg', 'ZeroOrderExclusionDeg', 'AbsoluteReflectancePercent'});
    group = struct('name', name, 'wavelengthNm', wavelength, 'centersDeg', centers, 'metrics', metrics);
end
