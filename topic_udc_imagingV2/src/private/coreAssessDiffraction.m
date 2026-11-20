function quality = coreAssessDiffraction(optics)
%coreAssessDiffraction 比较同一入瞳的遮挡与无遮挡成像。
    options = struct('tailRadiiUm',[10,20,50], 'mtfBandCyclesPerMm',[10,100], ...
        'hdrPointSignal',1000,'hdrDisplayHalfWidthUm',Inf);
    if isfield(optics.config,'quality')
        configured = optics.config.quality;
        for name = string(fieldnames(configured)).'
            assert(isfield(options,name),'UDCImagingV2:QualityOption','未知评价参数：%s',name);
            options.(name) = configured.(name);
        end
    end
    validateattributes(options.tailRadiiUm,{'numeric'},{'vector','real','finite','positive','increasing'});
    validateattributes(options.mtfBandCyclesPerMm,{'numeric'}, ...
        {'vector','numel',2,'real','finite','nonnegative','increasing'});
    validateattributes(options.hdrPointSignal,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(options.hdrDisplayHalfWidthUm,{'numeric'},{'scalar','real','nonnan','positive'});

    nChannels = size(optics.kernel,3);
    channels = (1:nChannels).';
    clearCaptured = reshape(sum(sum(optics.clearKernel,1),2),[],1);
    captured = optics.channelCapturedThroughput(:);
    quality.options = options;
    quality.channelMetrics = table(channels, optics.channelThroughput(:),captured,clearCaptured, ...
        divideWhenPositive(captured,clearCaptured),divideWhenPositive(captured,optics.channelThroughput(:)), ...
        'VariableNames',{'Channel','Throughput','CapturedThroughput','ClearCapturedThroughput', ...
        'RelativeToClear','FractionOfTransmittedCaptured'});

    clearPsf = optics.clearKernel./reshape(clearCaptured,1,1,[]);
    clearPsf(~isfinite(clearPsf)) = 0;
    otf = fft2(ifftshift(ifftshift(clearPsf,1),2));
    quality.clearPsf = clearPsf;
    quality.clearMtf = abs(fftshift(fftshift(otf,1),2));
    quality.tail = measureTail(optics,options.tailRadiiUm,captured,clearCaptured);
    [quality.mtfProfiles,quality.mtfSummary] = measureMtf(optics,quality.clearMtf,options.mtfBandCyclesPerMm);
    % 居中单位像素脉冲的响应就是能量核，直接缩放可保留线性动态范围。
    quality.hdrPoint = struct('clear',options.hdrPointSignal*optics.clearKernel, ...
        'withMask',options.hdrPointSignal*optics.kernel);
end

function metrics = measureTail(optics, radii, captured, clearCaptured)
    coordinateUm = optics.coordinate*1e6;
    distance = hypot(coordinateUm(:),coordinateUm);
    nChannels = numel(captured);
    tail = zeros(numel(radii),nChannels);
    clearTail = tail;
    for j = 1:numel(radii)
        outside = distance>radii(j);
        tail(j,:) = reshape(sum(sum(optics.kernel.*outside,1),2),1,[]);
        clearTail(j,:) = reshape(sum(sum(optics.clearKernel.*outside,1),2),1,[]);
    end
    fraction = divideWhenPositive(tail,repmat(captured.',numel(radii),1));
    clearFraction = divideWhenPositive(clearTail,repmat(clearCaptured.',numel(radii),1));
    % 圆形 ROI 超出窗口时，圆外份额不能作为该半径的有效对比指标。
    valid = radii(:)<=max(abs(coordinateUm));
    fraction(~valid,:) = NaN;
    clearFraction(~valid,:) = NaN;
    tail(~valid,:) = NaN;
    clearTail(~valid,:) = NaN;
    metrics = table(repelem((1:nChannels).',numel(radii),1),repmat(radii(:),nChannels,1), ...
        repmat(valid,nChannels,1),fraction(:),clearFraction(:),tail(:),clearTail(:), ...
        'VariableNames',{'Channel','RadiusUm','RadiusFitsWindow','TailFractionWithinWindow', ...
        'ClearTailFractionWithinWindow','TailThroughput','ClearTailThroughput'});
end

function [profiles, summary] = measureMtf(optics,clearMtf,band)
    axisFrequency = optics.spatialFrequency;
    frequency = axisFrequency(axisFrequency>=0);
    % 两端补到目标频段的精确频率，频段超出采样范围时仅保留曲线。
    covered = band(2)<=frequency(end);
    if covered, frequency = unique([frequency(:);band(:)]).'; end
    directions = ["Horizontal","Vertical","Diagonal45","Diagonal135"];
    angles = [0,90,45,135];
    nChannels = size(optics.mtf,3);
    profileRows = cell(nChannels,numel(angles));
    summaryRows = profileRows;
    for ch = 1:nChannels
        for j = 1:numel(angles)
            x = frequency*cosd(angles(j));
            y = frequency*sind(angles(j));
            values = interp2(axisFrequency,axisFrequency,optics.mtf(:,:,ch),x,y,'linear');
            reference = interp2(axisFrequency,axisFrequency,clearMtf(:,:,ch),x,y,'linear');
            profileRows{ch,j} = table(repmat(ch,numel(frequency),1),repmat(directions(j),numel(frequency),1), ...
                frequency(:),values(:),reference(:), ...
                'VariableNames',{'Channel','Direction','CyclesPerMm','MTF','ClearMTF'});
            minMtf = NaN;
            minClearMtf = NaN;
            if covered
                inBand = frequency>=band(1) & frequency<=band(2);
                minClearMtf = min(reference(inBand));
                if values(1)>0, minMtf = min(values(inBand)); end
            end
            summaryRows{ch,j} = table(ch,directions(j),firstHalfCrossing(frequency,values), ...
                firstHalfCrossing(frequency,reference),band(1),band(2),covered,minMtf,minClearMtf, ...
                'VariableNames',{'Channel','Direction','MTF50CyclesPerMm','ClearMTF50CyclesPerMm', ...
                'BandStartCyclesPerMm','BandEndCyclesPerMm','BandFitsAxis','BandMinMTF','ClearBandMinMTF'});
        end
    end
    profiles = vertcat(profileRows{:});
    summary = vertcat(summaryRows{:});
end

function crossing = firstHalfCrossing(frequency,values)
    crossing = NaN;
    idx = find(values<=0.5,1,'first');
    if ~isempty(idx) && idx>1
        fraction = (0.5-values(idx-1))/(values(idx)-values(idx-1));
        crossing = frequency(idx-1)+fraction*(frequency(idx)-frequency(idx-1));
    end
end

function quotient = divideWhenPositive(numerator,denominator)
    quotient = NaN(size(numerator));
    positive = denominator>0;
    quotient(positive) = numerator(positive)./denominator(positive);
end
