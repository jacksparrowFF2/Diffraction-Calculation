function spectrum = loadCoeSpectrum(wavelengthNm, dataDir)
%loadCoeSpectrum 读取 CIE 官方 D65 / 1931 2° 数据，生成光谱与色度积分权重。
% 仅对输入的波段积分；梯形积分包含端点半权重，避免遗漏或重复乘波长步长。
    arguments
        wavelengthNm (1,:) double {mustBeReal, mustBeFinite, mustBeInRange(wavelengthNm, 380, 780)}
        dataDir (1,1) string = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'Input', 'Colorimetry')
    end
    assert(numel(wavelengthNm) >= 2 && all(diff(wavelengthNm) > 0), ...
        'COEReflection:SpectrumGrid', '光谱计算至少需要两个严格递增的波长。');
    d65 = readmatrix(fullfile(dataDir, 'CIE_std_illum_D65.csv'));
    cmf = readmatrix(fullfile(dataDir, 'CIE_xyz_1931_2deg.csv'));
    wavelengths = wavelengthNm(:);
    spd = interp1(d65(:, 1), d65(:, 2), wavelengths, 'linear');
    xyz = interp1(cmf(:, 1), cmf(:, 2:4), wavelengths, 'linear');
    steps = diff(wavelengths);
    quadrature = [steps(1); steps(1:end-1) + steps(2:end); steps(end)]/2;
    rawWeights = spd.*quadrature;
    spectrum.wavelengthNm = wavelengths;
    spectrum.d65 = spd;
    spectrum.cmfXYZ = xyz;
    spectrum.quadratureNm = quadrature;
    spectrum.energyWeights = rawWeights/sum(rawWeights);
    % 理想单位反射体的 Y=1；所有 XYZ 使用同一常数，不能分别归一化三个通道。
    spectrum.xyzWeights = rawWeights.*xyz/sum(rawWeights.*xyz(:, 2));
    spectrum.photopicWeights = spectrum.xyzWeights(:, 2);
    spectrum.whiteXYZ = sum(spectrum.xyzWeights, 1);
    spectrum.whiteXY = spectrum.whiteXYZ(1:2)/sum(spectrum.whiteXYZ);
    spectrum.source = "CIE D65 + CIE 1931 2 degree; 380-780 nm band";
end
