function coefficient = reflectCoeFilmStack(wavelengths, incidentN, filmNK, thickness, bareReflection, referenceHeight)
%reflectCoeFilmStack 法向相干膜层递推；输入膜层从入射侧到下方排列，N=n+i*k。
% bareReflection 指同一底层直接接触入射介质时的复反射系数。
% r=-1 为理想电反射底；r=0 为与入射介质阻抗匹配的吸收终端（进入的能量不返回）。
% CF/BM 与该终端的界面仍按阻抗变换计算；不是把覆盖膜后的局部反射强制置零。
    arguments
        wavelengths (1,:) double {mustBeReal, mustBeFinite, mustBePositive}
        incidentN (1,1) double {mustBeReal, mustBeFinite, mustBePositive}
        filmNK double {mustBeFinite}
        thickness (1,:) double {mustBeReal, mustBeFinite, mustBeNonnegative}
        bareReflection double {mustBeFinite}
        referenceHeight (1,1) double {mustBeReal, mustBeFinite, mustBeNonnegative}
    end
    nLambda = numel(wavelengths);
    assert(size(filmNK, 1) == numel(thickness) && ...
        (isempty(filmNK) || size(filmNK, 2) == 1 || size(filmNK, 2) == nLambda), ...
        'COEReflection:FilmSize', '膜层 n/k 的行数应等于厚度数，列数为 1 或波长数。');
    assert(all(real(filmNK) > 0, 'all') && all(imag(filmNK) >= 0, 'all'), ...
        'COEReflection:PassiveNK', '膜层要求 n>0、k>=0。');
    assert(isscalar(bareReflection) || numel(bareReflection) == nLambda, ...
        'COEReflection:BareSize', '裸底反射系数必须为标量或逐波长向量。');
    assert(all(abs(bareReflection(:)) <= 1 + 1e-12), ...
        'COEReflection:PassiveBase', '裸底反射系数的模不得大于 1。');
    assert(referenceHeight + 1e-15 >= sum(thickness), ...
        'COEReflection:ReferencePlane', '共同参考平面不得低于膜堆顶部。');
    active = thickness > 0;
    filmNK = filmNK(active, :);
    thickness = thickness(active);
    if size(filmNK, 2) == 1
        filmNK = repmat(filmNK, 1, nLambda);
    end
    bare = reshape(bareReflection, 1, []).*ones(1, nLambda);
    coefficient = bare;
    if ~isempty(thickness)
        lastN = filmNK(end, :);
        % 将固定底层阻抗从入射介质变换到最下面的膜内；此写法对 r=-1 也无奇点。
        coefficient = ((lastN - incidentN) + (lastN + incidentN).*bare) ...
            ./((lastN + incidentN) + (lastN - incidentN).*bare);
        for layer = numel(thickness):-1:1
            previousN = incidentN;
            if layer > 1
                previousN = filmNK(layer-1, :);
            end
            interface = (previousN - filmNK(layer, :))./(previousN + filmNK(layer, :));
            roundTrip = coefficient.*exp(4i*pi*filmNK(layer, :)*thickness(layer)./wavelengths);
            coefficient = (interface + roundTrip)./(1 + interface.*roundTrip);
        end
    end
    % 各柱体的顶面高度不同，补足到同一空气参考平面的往返传播，不能分别丢弃这部分相位。
    airGap = max(0, referenceHeight - sum(thickness));
    coefficient = coefficient.*exp(4i*pi*incidentN*airGap./wavelengths);
    assert(all(isfinite(coefficient)) && all(abs(coefficient) <= 1 + 1e-9), ...
        'COEReflection:FilmPassivity', '膜层计算出现非有限值或无源反射率超过 1。');
end
