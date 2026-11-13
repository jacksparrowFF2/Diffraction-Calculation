function spectrum = loadCoeSpectrum(wavelengthNm, dataDir, spectrumFile, discreteLine)
%loadCoeSpectrum 读取 CIE 官方 D65 / 1931 2° 数据，生成光谱与色度积分权重。
% 仅对输入的波段积分；梯形积分包含端点半权重，避免遗漏或重复乘波长步长。
    arguments
        wavelengthNm (1,:) double {mustBeReal, mustBeFinite, mustBeInRange(wavelengthNm, 380, 780)}
        dataDir (1,1) string = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'Input', 'Colorimetry')
        spectrumFile (1,1) string = ""
        discreteLine (1,1) logical = false
    end
    assert((discreteLine && numel(wavelengthNm)==1) || ...
        (~discreteLine && numel(wavelengthNm)>=2 && all(diff(wavelengthNm)>0)), ...
        'COEReflection:SpectrumGrid','连续光谱需要至少两个递增波长；单色模式需要一个波长。');
    d65 = readmatrix(fullfile(dataDir, 'CIE_std_illum_D65.csv'));
    cmf = readmatrix(fullfile(dataDir, 'CIE_xyz_1931_2deg.csv'));
    wavelengths = wavelengthNm(:);
    spd = interp1(d65(:, 1), d65(:, 2), wavelengths, 'linear');
    source="CIE D65";
    if discreteLine
        spd=1; source="单色谱线";
    elseif strlength(spectrumFile)>0
        measured=readmatrix(spectrumFile);
        assert(size(measured,2)==2 && size(measured,1)>=2 && all(isfinite(measured),'all') ...
            && all(diff(measured(:,1))>0) && all(measured(:,2)>=0), ...
            'COEReflection:SourceSpectrum','光源CSV应为递增波长nm与非负相对光谱功率两列。');
        assert(measured(1,1)<=min(wavelengths) && measured(end,1)>=max(wavelengths), ...
            'COEReflection:SourceSpectrumRange','光源CSV必须覆盖选定的 wavelengthRangeNm，不进行外推。');
        spd=interp1(measured(:,1),measured(:,2),wavelengths,'linear');
        source="自定义光源";
    end
    xyz = interp1(cmf(:, 1), cmf(:, 2:4), wavelengths, 'linear');
    if discreteLine
        quadrature=1; % 单条谱线按总谱线功率归一化，不解释为每nm连续SPD
    else
        steps = diff(wavelengths);
        quadrature = [steps(1); steps(1:end-1) + steps(2:end); steps(end)]/2;
    end
    rawWeights = spd.*quadrature;
    assert(sum(rawWeights)>0 && sum(rawWeights.*xyz(:,2))>0, ...
        'COEReflection:SourceSpectrumPower','所选波段的光源功率和明视觉加权功率必须大于0。');
    spectrum.wavelengthNm = wavelengths;
    spectrum.d65 = spd;
    spectrum.incidentSPD=spd; % d65保留历史字段名；自定义/单色时以本字段和source辨别
    spectrum.discreteLine=discreteLine;
    spectrum.cmfXYZ = xyz;
    spectrum.quadratureNm = quadrature;
    spectrum.energyWeights = rawWeights/sum(rawWeights);
    % 理想单位反射体的 Y=1；所有 XYZ 使用同一常数，不能分别归一化三个通道。
    spectrum.xyzWeights = rawWeights.*xyz/sum(rawWeights.*xyz(:, 2));
    spectrum.photopicWeights = spectrum.xyzWeights(:, 2);
    spectrum.whiteXYZ = sum(spectrum.xyzWeights, 1);
    spectrum.whiteXY = spectrum.whiteXYZ(1:2)/sum(spectrum.whiteXYZ);
    spectrum.source = source + " + CIE 1931 2 degree";
end
