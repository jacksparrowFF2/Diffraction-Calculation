function lookup = coreTransmitLayerStack(wavelengths, stack, patterns)
%coreTransmitLayerStack 求各区域法向复透射，以入射功率归一化振幅。
    validateattributes(stack.incidentN, {'numeric'}, {'scalar','real','finite','positive'});
    validateattributes(stack.exitN, {'numeric'}, {'scalar','real','finite','positive'});
    names = string(stack.layerNames);
    direct = false(1,numel(names));
    directAmplitude = zeros(numel(names),numel(wavelengths));
    nk = zeros(numel(names),numel(wavelengths));
    thickness = zeros(1,numel(names));
    for j = 1:numel(names)
        layer = stack.layers.(names(j));
        mustBeMember(string(layer.method), ["transmittance","nk"]);
        direct(j) = string(layer.method) == "transmittance";
        if direct(j)
            validateattributes(layer.transmittance, {'numeric'}, {'vector','real','finite','>=',0,'<=',1});
            assert(isscalar(layer.transmittance) || numel(layer.transmittance)==numel(wavelengths), ...
                'UDCImagingV2:TransmittanceSize', '透过率应为标量或与波长一一对应。');
            directAmplitude(j,:) = sqrt(reshape(layer.transmittance,1,[]));
        else
            validateattributes(layer.thickness, {'numeric'}, {'scalar','real','finite','nonnegative'});
            validateattributes(layer.n, {'numeric'}, {'vector','real','finite','positive'});
            validateattributes(layer.k, {'numeric'}, {'vector','real','finite','nonnegative'});
            assert(isscalar(layer.n) || numel(layer.n)==numel(wavelengths), ...
                'UDCImagingV2:IndexSize', 'n 应为标量或与波长一一对应。');
            assert(isscalar(layer.k) || numel(layer.k)==numel(wavelengths), ...
                'UDCImagingV2:IndexSize', 'k 应为标量或与波长一一对应。');
            nk(j,:) = reshape(layer.n,1,[]) + 1i*reshape(layer.k,1,[]);
            thickness(j) = layer.thickness;
        end
    end
    referenceHeight = 0;
    for level = unique(stack.layerLevels)
        referenceHeight = referenceHeight + max(thickness(stack.layerLevels==level));
    end

    lookup = complex(zeros(size(patterns,1),numel(wavelengths)));
    for region = 1:size(patterns,1)
        present = find(patterns(region,:));
        if isempty(present)
            continue
        end
        effective = present(direct(present));
        if ~isempty(effective)
            assert(isscalar(present), 'UDCImagingV2:DirectOverlap', ...
                '直接透过率代表区域的整体透射，白区不可与其他透射层重叠。');
            lookup(region,:) = directAmplitude(effective,:);
            continue
        end
        % 图层从下到上存储，从出射界面向入射侧递推。
        films = find(patterns(region,:) & thickness>0);
        indices = [stack.exitN*ones(size(wavelengths)); nk(films,:); stack.incidentN*ones(size(wavelengths))];
        r = zeros(size(wavelengths));
        t = ones(size(wavelengths));
        p = t;
        for j = 2:size(indices,1)
            lower = indices(j-1,:);
            upper = indices(j,:);
            interface = (upper-lower)./(upper+lower);
            denominator = 1 + interface.*r.*p.^2;
            t = (1+interface).*p.*t./denominator;
            r = (interface+r.*p.^2)./denominator;
            if j <= numel(films)+1
                p = exp(2i*pi*upper*thickness(films(j-1))./wavelengths);
            end
        end
        gap = max(0,referenceHeight-sum(thickness(films)));
        lookup(region,:) = sqrt(stack.exitN/stack.incidentN)*t ...
            .*exp(2i*pi*stack.incidentN*gap./wavelengths);
    end
    assert(all(isfinite(lookup),'all') && all(abs(lookup)<=1+1e-9,'all'), ...
        'UDCImagingV2:FilmPassivity', '膜堆透射系数须有限且满足无源功率约束。');
end
