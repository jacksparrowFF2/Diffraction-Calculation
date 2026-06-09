%% 多形状相位微岛 DOE 优化
% 本脚本用于搜索“非零级衍射峰最小 + 0 级主斑保持 + 少衍射环 + haze 风险受限”的高折微岛设计。
% 当前支持 PatternMode:
%   circle：多半径圆形相位微岛 baseline；
%   random_ellipse：多长轴、随机旋角椭圆相位微岛；
%   blue_noise_circle / blue_noise_ellipse：best-candidate 蓝噪声中心分布。
% 优化逻辑和可视化脚本分离：本脚本负责筛选候选设计，
% CustomApertureWithPhase_fraunhofer_fft.m 负责对某一个设计做详细出图。

clear, clc, close all

%% DOE 运行模式
% quick：用于快速调通流程；full：用于正式粗扫，seed 使用 1:10。
doePreset = "quick";

switch doePreset
    case "quick"
        nn = 180;
        seedList = 1:2;
        radiusSetList_um = {
            [1.0, 1.5]
            [0.8, 1.2, 1.8]
            [0.8, 1.2, 1.6, 2.0]
        };        
        patternModeList = ["circle", "random_ellipse", "blue_noise_circle", "blue_noise_ellipse"];
        ellipseAspectRatioList = [1.5];
        phaseHeightList_nm = [140, 160, 180];
        fillFactorList = [0.10, 0.15];
        minDistanceFactorList = [0.8, 1.0, 1.2];
    case "full"
        nn = 320;
        seedList = 1:10;
        radiusSetList_um = {
            [1.0, 1.5]
            [0.8, 1.2, 1.8]
            [0.8, 1.2, 1.6, 2.0]
            [0.6, 1.0, 1.4, 1.8]
            [0.8, 1.0, 1.4, 1.8, 2.2]
        };
        patternModeList = ["circle", "random_ellipse", "blue_noise_circle", "blue_noise_ellipse"];
        ellipseAspectRatioList = [1.5, 2.0, 2.5];
        phaseHeightList_nm = [60, 100, 140, 180, 220];
        fillFactorList = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30];
        minDistanceFactorList = [0.8, 1.0, 1.2, 1.5, 2.0];
    otherwise
        error('未知 DOE 预设：%s', doePreset);
end

%% 基础光学与几何参数
% 固定高度、折射率差和双程相位，与当前可视化脚本保持一致。
lambdalist = [620e-9, 550e-9, 450e-9];
z = 0.2;
objectunit = 1e-6;
Imageunit = 1e-3;
Pitch = 50 * objectunit;
r = 10 * objectunit;
Gratingflag = 1;

Nx_period = 2;
Ny_period = 2;
Fx = 20 * Imageunit;
Fy = 20 * Imageunit;
Angle = 2.5;
FocusR = tand(Angle) * z;

n_high = 1.75;
n_low = 1.55;
delta_n = n_high - n_low;
phaseHeightList = phaseHeightList_nm * 1e-9;
patternConfigList = buildPatternConfigList(radiusSetList_um, patternModeList, ellipseAspectRatioList);

% 目标函数权重。0 级主斑、近角衍射环和 haze 风险都作为强约束处理。
topK = 12;
hazeLimitMultiplier = 1.10;
nearRingLimitMultiplier = 1.03;
ringVisibilityLimitMultiplier = 1.03;
zeroOrderRadiusFactor = 0.45;
zeroOrderRetentionTarget = 0.98;
radialBinCount = 90;
weightPeakBackground = 0.15;
weightZeroOrderPenalty = 10.0;
weightNearRingPenalty = 8.0;
weightRingVisibilityPenalty = 8.0;
weightHazePenalty = 10.0;

%% 物面和像面网格
x_cell = linspace(-Pitch/2, Pitch/2, nn);
y_cell = linspace(-Pitch/2, Pitch/2, nn);
[Xc, Yc] = meshgrid(x_cell, y_cell);

xmin = -Nx_period * Pitch / 2;
xmax = -xmin;
ymin = -Ny_period * Pitch / 2;
ymax = -ymin;

ImageMaskRatioX = Fx / xmax;
ImageMaskRatioY = Fy / ymax;
Xmin = xmin * ImageMaskRatioX;
Xmax = -Xmin;
Ymin = ymin * ImageMaskRatioY;
Ymax = -Ymin;

Tnn = Nx_period * nn;
xVector = linspace(xmin, xmax, Tnn);
yVector = linspace(ymin, ymax, Tnn);
[x, y] = meshgrid(xVector, yVector);
XVector = linspace(Xmin, Xmax, Tnn);
YVector = linspace(Ymin, Ymax, Tnn);
[X, Y] = meshgrid(XVector, YVector);

%% 构建基础周期孔径
if Gratingflag == 1
    cellunit = 2;
    SinglePattern1 = CustomPattern(30, Xc, Yc, r, r, 0, 0);
    SinglePattern2 = CustomPattern(0, Xc, Yc, r, r, 0, 0);
    SinglePattern3 = CustomPattern(0, Xc, Yc, r, r, 0, 0);
    SinglePattern4 = CustomPattern(-30, Xc, Yc, r, r, 0, 0);
    apertureCell = [SinglePattern1 SinglePattern2; SinglePattern3 SinglePattern4];
    BaseMask = repmat(apertureCell, Ny_period/cellunit, Nx_period/cellunit);
else
    error('当前优化脚本先针对 Gratingflag == 1 的四孔径晶胞。');
end

%% 预计算无相位层基准指标
fprintf('开始预计算无相位层基准远场...\n');
baseline = repmat(struct(), 1, numel(lambdalist));
for c = 1:numel(lambdalist)
    lambda = lambdalist(c);
    k = 2 * pi / lambda;
    [~, U2_NoPhase] = fraunhofer_fft_SPH(BaseMask, xmin, xmax, ymin, ymax, Tnn, Tnn, ...
        lambda, z, Xmin, Xmax, Ymin, Ymax, Tnn, Tnn, k, 0);
    I0 = abs(U2_NoPhase).^2;
    zeroOrderRadius = zeroOrderRadiusFactor * lambda * z / Pitch;
    zeroOrderMask = X.^2 + Y.^2 <= zeroOrderRadius^2;
    baseline(c).I = I0;
    baseline(c).zeroOrderEnergy = sum(I0(zeroOrderMask), 'all');
    baseline(c).topKMean = meanTopKLocalPeaks(I0, topK, zeroOrderMask);
    baseline(c).peakBackground = peakBackgroundRatio(I0, topK, zeroOrderMask);
    baseline(c).nearRingEnergy = nearRingEnergyRatio(I0, X, Y, zeroOrderRadius, FocusR);
    baseline(c).ringVisibility = ringVisibilityByRadius(I0, X, Y, zeroOrderRadius, FocusR, radialBinCount);
    baseline(c).hazeRisk = hazeRiskByAngle(I0, X, Y, FocusR);
end
baselineHazeMean = mean([baseline.hazeRisk]);
hazeRiskLimit = baselineHazeMean * hazeLimitMultiplier;
baselineNearRingMean = mean([baseline.nearRingEnergy]);
nearRingEnergyLimit = baselineNearRingMean * nearRingLimitMultiplier;
baselineRingVisibilityMean = mean([baseline.ringVisibility]);
ringVisibilityLimit = baselineRingVisibilityMean * ringVisibilityLimitMultiplier;
fprintf('无相位层平均 2.5° 外能量比例：%.4f，haze 风险上限：%.4f\n', ...
    baselineHazeMean, hazeRiskLimit);
fprintf('无相位层近角环能量比例：%.4f，上限：%.4f\n', ...
    baselineNearRingMean, nearRingEnergyLimit);
fprintf('无相位层径向环可见度：%.4f，上限：%.4f\n', ...
    baselineRingVisibilityMean, ringVisibilityLimit);

%% DOE 粗扫
totalDesigns = numel(patternConfigList) * numel(phaseHeightList) * numel(fillFactorList) * ...
    numel(minDistanceFactorList) * numel(seedList);
fprintf('开始 DOE：%d 个设计，模式：%s\n', totalDesigns, doePreset);

results = table();
designIndex = 0;
bestObjective = inf;
bestDesign = struct();

for patternConfigId = 1:numel(patternConfigList)
    patternConfig = patternConfigList(patternConfigId);
    sizeList = patternConfig.SizeSet_um * 1e-6;
    distanceRadius = patternConfig.DistanceRadius_um * 1e-6;

    for phaseHeight = phaseHeightList
        for fillFactor = fillFactorList
            for minDistanceFactor = minDistanceFactorList
                minCenterDistance = minDistanceFactor * 2 * distanceRadius;

                for seed = seedList
                    designIndex = designIndex + 1;
                    fprintf('[%d/%d] mode=%s, sizeSet=%d, aspect=%.2f, height=%.0f nm, fill=%.2f, minDistFactor=%.2f, seed=%d\n', ...
                        designIndex, totalDesigns, patternConfig.PatternMode, patternConfig.SizeSetId, ...
                        patternConfig.AspectRatio, phaseHeight * 1e9, fillFactor, minDistanceFactor, seed);

                    [h_map, phaseCenters, actualFillFactor] = generatePoissonDiskPhaseHeightMap( ...
                        x, y, xmin, xmax, ymin, ymax, phaseHeight, patternConfig.PatternMode, ...
                        sizeList, patternConfig.AspectRatio, minCenterDistance, fillFactor, seed);
                    fillAchievementRatio = actualFillFactor / max(fillFactor, eps);
                    fillShortfall = max(0, fillFactor - actualFillFactor);
                    effectiveFillBin = round(actualFillFactor / 0.01) * 0.01;

                    metric = evaluatePhaseDesign(BaseMask, h_map, baseline, lambdalist, ...
                        delta_n, xmin, xmax, ymin, ymax, Tnn, z, Xmin, Xmax, Ymin, Ymax, ...
                        X, Y, FocusR, Pitch, topK, zeroOrderRadiusFactor, zeroOrderRetentionTarget, ...
                        radialBinCount, nearRingEnergyLimit, ringVisibilityLimit, hazeRiskLimit, ...
                        weightPeakBackground, weightZeroOrderPenalty, weightNearRingPenalty, ...
                        weightRingVisibilityPenalty, weightHazePenalty);

                    newRow = table( ...
                        designIndex, patternConfigId, patternConfig.SizeSetId, string(patternConfig.PatternMode), ...
                        string(mat2str(patternConfig.SizeSet_um)), patternConfig.AspectRatio, ...
                        phaseHeight * 1e9, fillFactor, minDistanceFactor, seed, size(phaseCenters, 1), ...
                        actualFillFactor, fillAchievementRatio, fillShortfall, effectiveFillBin, ...
                        metric.NonZeroTopKPeakRatio, metric.PeakBackgroundRatioPenalty, ...
                        metric.ZeroOrderRetention, metric.ZeroOrderLossPenalty, ...
                        metric.NearRingEnergyRatio, metric.NearRingEnergyPenalty, ...
                        metric.RingVisibility, metric.RingVisibilityRatio, metric.RingVisibilityPenalty, ...
                        metric.HazeRisk, metric.HazeRiskPenalty, metric.Objective, ...
                        'VariableNames', {'DesignIndex', 'PatternConfigId', 'SizeSetId', 'PatternMode', ...
                        'RadiusSet_um', 'AspectRatio', 'PhaseHeight_nm', ...
                        'TargetFillFactor', 'MinDistanceFactor', 'Seed', ...
                        'IslandCount', 'ActualFillFactor', 'FillAchievementRatio', ...
                        'FillShortfall', 'EffectiveFillBin', 'NonZeroTopKPeakRatio', ...
                        'PeakBackgroundRatioPenalty', 'ZeroOrderRetention', 'ZeroOrderLossPenalty', ...
                        'NearRingEnergyRatio', 'NearRingEnergyPenalty', ...
                        'RingVisibility', 'RingVisibilityRatio', 'RingVisibilityPenalty', ...
                        'HazeRisk', 'HazeRiskPenalty', 'Objective'});
                    results = [results; newRow]; %#ok<AGROW>

                    if metric.Objective < bestObjective
                        bestObjective = metric.Objective;
                        bestDesign.h_map = h_map;
                        bestDesign.centers = phaseCenters;
                        bestDesign.patternMode = patternConfig.PatternMode;
                        bestDesign.sizeSet_um = patternConfig.SizeSet_um;
                        bestDesign.aspectRatio = patternConfig.AspectRatio;
                        bestDesign.radiusSet_um = patternConfig.SizeSet_um;
                        bestDesign.phaseHeight = phaseHeight;
                        bestDesign.fillFactor = fillFactor;
                        bestDesign.minDistanceFactor = minDistanceFactor;
                        bestDesign.seed = seed;
                        bestDesign.metric = metric;
                        bestDesign.actualFillFactor = actualFillFactor;
                    end
                end
            end
        end
    end
end

results = sortrows(results, 'Objective', 'ascend');
disp(results(1:min(10, height(results)), :));

%% 保存结果
outputDir = fullfile('..', 'ImageForShow');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

csvPath = fullfile(outputDir, 'multi_radius_phase_doe_results.csv');
matPath = fullfile(outputDir, 'multi_radius_phase_doe_results.mat');
writetable(results, csvPath);
save(matPath, 'results', 'bestDesign', 'baselineHazeMean', 'hazeRiskLimit', ...
    'baselineNearRingMean', 'nearRingEnergyLimit', ...
    'baselineRingVisibilityMean', 'ringVisibilityLimit', 'doePreset');

figure;
imagesc(x(1,:) * 1e6, y(:,1) * 1e6, bestDesign.h_map * 1e9);
axis image;
set(gca, 'YDir', 'normal');
xlabel('x / um');
ylabel('y / um');
title(sprintf('Best Multi-radius Phase Height Map, Obj=%.4f', bestDesign.metric.Objective));
colormap gray;
colorbar;
exportgraphics(gcf, fullfile(outputDir, 'best_multi_radius_phase_height_map.png'), 'Resolution', 200);

fprintf('\n最佳设计：\n');
fprintf('patternMode = %s\n', bestDesign.patternMode);
fprintf('sizeSet_um = %s, aspectRatio = %.2f\n', mat2str(bestDesign.sizeSet_um), bestDesign.aspectRatio);
fprintf('phaseHeight = %.0f nm\n', bestDesign.phaseHeight * 1e9);
fprintf('fillFactor = %.2f, minDistanceFactor = %.2f, seed = %d\n', ...
    bestDesign.fillFactor, bestDesign.minDistanceFactor, bestDesign.seed);
fprintf('actualFillFactor = %.4f\n', bestDesign.actualFillFactor);
fprintf('NonZeroTopKPeakRatio = %.4f\n', bestDesign.metric.NonZeroTopKPeakRatio);
fprintf('PeakBackgroundRatioPenalty = %.4f\n', bestDesign.metric.PeakBackgroundRatioPenalty);
fprintf('ZeroOrderRetention = %.4f\n', bestDesign.metric.ZeroOrderRetention);
fprintf('NearRingEnergyRatio = %.4f\n', bestDesign.metric.NearRingEnergyRatio);
fprintf('RingVisibilityRatio = %.4f\n', bestDesign.metric.RingVisibilityRatio);
fprintf('HazeRisk = %.4f, HazeLimit = %.4f\n', bestDesign.metric.HazeRisk, hazeRiskLimit);
fprintf('Objective = %.4f\n', bestDesign.metric.Objective);
fprintf('结果已保存：%s\n', csvPath);

%% 局部函数
function patternConfigList = buildPatternConfigList(sizeSetList_um, patternModeList, ellipseAspectRatioList)
    %BUILDPATTERNCONFIGLIST 展开圆形和随机旋角椭圆的 DOE pattern 配置。
    patternConfigList = struct('PatternMode', {}, 'SizeSetId', {}, 'SizeSet_um', {}, ...
        'AspectRatio', {}, 'DistanceRadius_um', {});

    for mode = patternModeList
        switch string(mode)
            case {"circle", "blue_noise_circle"}
                for sizeSetId = 1:numel(sizeSetList_um)
                    patternConfigList(end + 1) = struct( ... %#ok<AGROW>
                        'PatternMode', string(mode), ...
                        'SizeSetId', sizeSetId, ...
                        'SizeSet_um', sizeSetList_um{sizeSetId}, ...
                        'AspectRatio', 1.0, ...
                        'DistanceRadius_um', max(sizeSetList_um{sizeSetId}));
                end

            case {"random_ellipse", "blue_noise_ellipse"}
                for sizeSetId = 1:numel(sizeSetList_um)
                    for aspectRatio = ellipseAspectRatioList
                        patternConfigList(end + 1) = struct( ... %#ok<AGROW>
                            'PatternMode', string(mode), ...
                            'SizeSetId', sizeSetId, ...
                            'SizeSet_um', sizeSetList_um{sizeSetId}, ...
                            'AspectRatio', aspectRatio, ...
                            'DistanceRadius_um', max(sizeSetList_um{sizeSetId}) / sqrt(aspectRatio));
                    end
                end

            otherwise
                error('未知 PatternMode：%s', string(mode));
        end
    end
end

function metric = evaluatePhaseDesign(BaseMask, h_map, baseline, lambdalist, ...
    delta_n, xmin, xmax, ymin, ymax, Tnn, z, Xmin, Xmax, Ymin, Ymax, ...
    X, Y, FocusR, Pitch, topK, zeroOrderRadiusFactor, zeroOrderRetentionTarget, ...
    radialBinCount, nearRingEnergyLimit, ringVisibilityLimit, hazeRiskLimit, ...
    weightPeakBackground, weightZeroOrderPenalty, weightNearRingPenalty, ...
    weightRingVisibilityPenalty, weightHazePenalty)
    %EVALUATEPHASEDESIGN 计算单个相位层设计的目标函数。

    topKRatioList = zeros(1, numel(lambdalist));
    peakBackgroundRatioList = zeros(1, numel(lambdalist));
    zeroOrderRetentionList = zeros(1, numel(lambdalist));
    nearRingEnergyList = zeros(1, numel(lambdalist));
    ringVisibilityList = zeros(1, numel(lambdalist));
    ringVisibilityRatioList = zeros(1, numel(lambdalist));
    hazeRiskList = zeros(1, numel(lambdalist));

    for c = 1:numel(lambdalist)
        lambda = lambdalist(c);
        k = 2 * pi / lambda;
        phi_double = 4 * pi / lambda * delta_n * h_map;
        ComplexMask = BaseMask .* exp(1j * phi_double);

        [~, U2_WithPhase] = fraunhofer_fft_SPH(ComplexMask, xmin, xmax, ymin, ymax, Tnn, Tnn, ...
            lambda, z, Xmin, Xmax, Ymin, Ymax, Tnn, Tnn, k, 0);
        I = abs(U2_WithPhase).^2;
        zeroOrderRadius = zeroOrderRadiusFactor * lambda * z / Pitch;
        zeroOrderMask = X.^2 + Y.^2 <= zeroOrderRadius^2;

        topKMean = meanTopKLocalPeaks(I, topK, zeroOrderMask);
        peakBackground = peakBackgroundRatio(I, topK, zeroOrderMask);
        zeroOrderEnergy = sum(I(zeroOrderMask), 'all');
        nearRingEnergy = nearRingEnergyRatio(I, X, Y, zeroOrderRadius, FocusR);
        ringVisibility = ringVisibilityByRadius(I, X, Y, zeroOrderRadius, FocusR, radialBinCount);
        hazeRisk = hazeRiskByAngle(I, X, Y, FocusR);

        topKRatioList(c) = topKMean / baseline(c).topKMean;
        peakBackgroundRatioList(c) = peakBackground / baseline(c).peakBackground;
        zeroOrderRetentionList(c) = zeroOrderEnergy / baseline(c).zeroOrderEnergy;
        nearRingEnergyList(c) = nearRingEnergy;
        ringVisibilityList(c) = ringVisibility;
        ringVisibilityRatioList(c) = ringVisibility / max(baseline(c).ringVisibility, eps);
        hazeRiskList(c) = hazeRisk;
    end

    metric.NonZeroTopKPeakRatio = mean(topKRatioList);
    metric.TopKPeakRatio = metric.NonZeroTopKPeakRatio;
    metric.PeakBackgroundRatioPenalty = mean(peakBackgroundRatioList);
    metric.ZeroOrderRetention = mean(zeroOrderRetentionList);
    metric.ZeroOrderLossPenalty = max(0, zeroOrderRetentionTarget - metric.ZeroOrderRetention)^2;
    metric.NearRingEnergyRatio = mean(nearRingEnergyList);
    metric.NearRingEnergyPenalty = max(0, metric.NearRingEnergyRatio - nearRingEnergyLimit)^2;
    metric.RingVisibility = mean(ringVisibilityList);
    metric.RingVisibilityRatio = mean(ringVisibilityRatioList);
    metric.RingVisibilityPenalty = max(0, metric.RingVisibility - ringVisibilityLimit)^2;
    metric.HazeRisk = mean(hazeRiskList);
    metric.HazeRiskPenalty = max(0, metric.HazeRisk - hazeRiskLimit)^2;
    metric.Objective = metric.NonZeroTopKPeakRatio + ...
        weightPeakBackground * metric.PeakBackgroundRatioPenalty + ...
        weightZeroOrderPenalty * metric.ZeroOrderLossPenalty + ...
        weightNearRingPenalty * metric.NearRingEnergyPenalty + ...
        weightRingVisibilityPenalty * metric.RingVisibilityPenalty + ...
        weightHazePenalty * metric.HazeRiskPenalty;
end

function [hMap, centers, actualFillFactor] = generatePoissonDiskPhaseHeightMap( ...
    xGrid, yGrid, xmin, xmax, ymin, ymax, islandHeight, patternMode, sizeList, ...
    aspectRatio, minCenterDistance, targetFillFactor, randomSeed)
    %GENERATEPOISSONDISKPHASEHEIGHTMAP 生成泊松盘圆形或随机旋角椭圆微岛高度图。
    % centers = [cx, cy, majorAxis, minorAxis, rotationAngleRad]。

    rng(randomSeed);
    xVector = xGrid(1, :);
    yVector = yGrid(:, 1);
    hMap = zeros(size(xGrid));

    majorAxisList = sizeList(:).';
    switch string(patternMode)
        case {"circle", "blue_noise_circle"}
            minorAxisList = majorAxisList;
            shapeAreaList = pi * majorAxisList.^2;
        case {"random_ellipse", "blue_noise_ellipse"}
            minorAxisList = majorAxisList ./ aspectRatio;
            shapeAreaList = pi * majorAxisList .* minorAxisList;
        otherwise
            error('未知 PatternMode：%s', string(patternMode));
    end

    maxMajorAxis = max(majorAxisList);
    validXmin = xmin + maxMajorAxis;
    validXmax = xmax - maxMajorAxis;
    validYmin = ymin + maxMajorAxis;
    validYmax = ymax - maxMajorAxis;
    if validXmin >= validXmax || validYmin >= validYmax
        error('微岛最大长半轴过大，已经超过当前物面尺寸。');
    end

    domainArea = (xmax - xmin) * (ymax - ymin);
    meanIslandArea = mean(shapeAreaList);
    targetIslandCount = ceil(targetFillFactor * domainArea / meanIslandArea);
    maxAttempts = max(20000, targetIslandCount * 1000);

    centers = zeros(targetIslandCount, 5);
    centerCount = 0;
    coveredArea = 0;
    targetArea = targetFillFactor * domainArea;
    attemptCount = 0;
    useBlueNoisePlacement = startsWith(string(patternMode), "blue_noise");
    blueNoiseCandidateCount = 24;
    blueNoiseSoftDistanceFactor = 0.65;

    while coveredArea < targetArea && attemptCount < maxAttempts
        attemptCount = attemptCount + 1;
        if useBlueNoisePlacement
            [candidateX, candidateY, sizeIndex, nearestDistanceSquared] = chooseBlueNoiseCandidate( ...
                centers, centerCount, validXmin, validXmax, validYmin, validYmax, ...
                majorAxisList, blueNoiseCandidateCount);
        else
            sizeIndex = randi(numel(majorAxisList));
            candidateX = validXmin + rand() * (validXmax - validXmin);
            candidateY = validYmin + rand() * (validYmax - validYmin);
            nearestDistanceSquared = inf;
        end
        candidateMajorAxis = majorAxisList(sizeIndex);
        candidateMinorAxis = minorAxisList(sizeIndex);
        candidateAngle = pi * rand();

        if centerCount == 0
            acceptCandidate = true;
        elseif useBlueNoisePlacement
            softMinDistance = blueNoiseSoftDistanceFactor * minCenterDistance;
            acceptCandidate = nearestDistanceSquared >= softMinDistance^2;
        else
            dx = centers(1:centerCount, 1) - candidateX;
            dy = centers(1:centerCount, 2) - candidateY;
            acceptCandidate = all(dx.^2 + dy.^2 >= minCenterDistance^2);
        end

        if acceptCandidate
            centerCount = centerCount + 1;
            if centerCount > size(centers, 1)
                centers(end + targetIslandCount, :) = 0; %#ok<AGROW>
            end
            centers(centerCount, :) = [candidateX, candidateY, candidateMajorAxis, candidateMinorAxis, candidateAngle];
            coveredArea = coveredArea + pi * candidateMajorAxis * candidateMinorAxis;
        end
    end

    centers = centers(1:centerCount, :);

    for idx = 1:centerCount
        cx = centers(idx, 1);
        cy = centers(idx, 2);
        majorAxis = centers(idx, 3);
        minorAxis = centers(idx, 4);
        rotationAngle = centers(idx, 5);

        xIndex = find(abs(xVector - cx) <= majorAxis);
        yIndex = find(abs(yVector - cy) <= majorAxis);
        [localX, localY] = meshgrid(xVector(xIndex), yVector(yIndex));
        dx = localX - cx;
        dy = localY - cy;
        localMajorCoord = dx * cos(rotationAngle) + dy * sin(rotationAngle);
        localMinorCoord = -dx * sin(rotationAngle) + dy * cos(rotationAngle);
        localShape = (localMajorCoord / majorAxis).^2 + (localMinorCoord / minorAxis).^2 <= 1;

        localHeight = hMap(yIndex, xIndex);
        localHeight(localShape) = islandHeight;
        hMap(yIndex, xIndex) = localHeight;
    end

    actualFillFactor = nnz(hMap > 0) / numel(hMap);
end

function [candidateX, candidateY, sizeIndex, nearestDistanceSquared] = chooseBlueNoiseCandidate( ...
    centers, centerCount, validXmin, validXmax, validYmin, validYmax, majorAxisList, candidateCount)
    %CHOOSEBLUENOISECANDIDATE Mitchell best-candidate 风格的蓝噪声中心选点。
    candidateXList = validXmin + rand(candidateCount, 1) * (validXmax - validXmin);
    candidateYList = validYmin + rand(candidateCount, 1) * (validYmax - validYmin);
    sizeIndexList = randi(numel(majorAxisList), candidateCount, 1);

    if centerCount == 0
        scoreList = min([candidateXList - validXmin, validXmax - candidateXList, ...
            candidateYList - validYmin, validYmax - candidateYList], [], 2).^2;
    else
        scoreList = zeros(candidateCount, 1);
        for idx = 1:candidateCount
            dx = centers(1:centerCount, 1) - candidateXList(idx);
            dy = centers(1:centerCount, 2) - candidateYList(idx);
            centerDistanceSquared = min(dx.^2 + dy.^2);
            boundaryDistance = min([candidateXList(idx) - validXmin, validXmax - candidateXList(idx), ...
                candidateYList(idx) - validYmin, validYmax - candidateYList(idx)]);
            scoreList(idx) = min(centerDistanceSquared, boundaryDistance^2);
        end
    end

    [nearestDistanceSquared, bestIndex] = max(scoreList);
    candidateX = candidateXList(bestIndex);
    candidateY = candidateYList(bestIndex);
    sizeIndex = sizeIndexList(bestIndex);
end

function value = meanTopKLocalPeaks(I, topK, excludeMask)
    %MEANTOPKLOCALPEAKS 提取局部峰并计算前 K 个峰的平均强度。
    if nargin < 3
        excludeMask = false(size(I));
    end
    localMaxMask = true(size(I));
    for rowShift = -1:1
        for colShift = -1:1
            if rowShift == 0 && colShift == 0
                continue
            end
            localMaxMask = localMaxMask & I >= shiftWithNegInf(I, rowShift, colShift);
        end
    end
    localMaxMask([1 end], :) = false;
    localMaxMask(:, [1 end]) = false;
    localMaxMask = localMaxMask & I > 0 & ~excludeMask;
    peakValues = sort(I(localMaxMask), 'descend');
    if isempty(peakValues)
        value = max(I(:));
        return
    end
    k = min(topK, numel(peakValues));
    value = mean(peakValues(1:k));
end

function shifted = shiftWithNegInf(I, rowShift, colShift)
    %SHIFTWITHNEGINF 平移矩阵，空出的边界填 -Inf，避免边界环绕。
    shifted = -inf(size(I));
    rowCount = size(I, 1);
    colCount = size(I, 2);

    sourceRows = max(1, 1 - rowShift):min(rowCount, rowCount - rowShift);
    sourceCols = max(1, 1 - colShift):min(colCount, colCount - colShift);
    targetRows = sourceRows + rowShift;
    targetCols = sourceCols + colShift;

    shifted(targetRows, targetCols) = I(sourceRows, sourceCols);
end

function ratio = peakBackgroundRatio(I, topK, excludeMask)
    %PEAKBACKGROUNDRATIO 用前 K 个局部峰均值除以背景中位数，估计彩色峰可见风险。
    if nargin < 3
        excludeMask = false(size(I));
    end
    peakValue = meanTopKLocalPeaks(I, topK, excludeMask);
    backgroundSamples = I(~excludeMask);
    threshold = prctile(backgroundSamples(:), 90);
    background = median(backgroundSamples(backgroundSamples <= threshold));
    ratio = peakValue / max(background, eps);
end

function ratio = nearRingEnergyRatio(I, X, Y, innerRadius, outerRadius)
    %NEARRINGENERGYRATIO 统计 0 级主斑外、haze 角度内的能量比例。
    R2 = X.^2 + Y.^2;
    ringMask = R2 > innerRadius^2 & R2 <= outerRadius^2;
    ratio = sum(I(ringMask), 'all') / max(sum(I, 'all'), eps);
end

function visibility = ringVisibilityByRadius(I, X, Y, innerRadius, outerRadius, binCount)
    %RINGVISIBILITYBYRADIUS 用径向平均曲线估计同心衍射环可见度。
    % 指标由局部径向峰的相对 prominence 和峰数量组成；数值越大，环状结构越明显。
    R = sqrt(X.^2 + Y.^2);
    ringMask = R > innerRadius & R <= outerRadius;
    if ~any(ringMask, 'all')
        visibility = 0;
        return
    end

    radiusSamples = R(ringMask);
    binIndex = floor((radiusSamples - innerRadius) / (outerRadius - innerRadius) * binCount) + 1;
    binIndex = min(max(binIndex, 1), binCount);
    intensitySamples = I(ringMask);
    validMask = isfinite(binIndex);
    binIndex = binIndex(validMask);
    intensitySamples = intensitySamples(validMask);

    radialSum = accumarray(binIndex(:), intensitySamples(:), [binCount, 1], @sum, 0);
    radialCount = accumarray(binIndex(:), 1, [binCount, 1], @sum, 0);
    radialMean = radialSum ./ max(radialCount, 1);
    radialMean = fillmissing(radialMean, 'linear', 'EndValues', 'nearest');
    smoothMean = movingAverageSame(radialMean, 7);
    residual = radialMean - smoothMean;

    peakMask = false(size(radialMean));
    peakMask(2:end-1) = residual(2:end-1) > 0 & ...
        radialMean(2:end-1) > radialMean(1:end-2) & ...
        radialMean(2:end-1) >= radialMean(3:end);

    positiveProminence = max(residual(peakMask), 0);
    meanLevel = mean(radialMean);
    peakProminenceScore = sum(positiveProminence) / max(meanLevel, eps);
    peakCountScore = 0.03 * nnz(peakMask);
    visibility = peakProminenceScore + peakCountScore;
end

function y = movingAverageSame(x, windowSize)
    %MOVINGAVERAGESAME 不依赖工具箱的同长度移动平均。
    windowSize = max(1, round(windowSize));
    kernel = ones(windowSize, 1) / windowSize;
    y = conv(x(:), kernel, 'same');
end

function hazeRisk = hazeRiskByAngle(I, X, Y, focusRadius)
    %HAZERISKBYANGLE 使用 2.5° 对应半径外的能量比例近似 haze 风险。
    angleMask = X.^2 + Y.^2 <= focusRadius^2;
    hazeRisk = sum(I(~angleMask), 'all') / sum(I, 'all');
end
