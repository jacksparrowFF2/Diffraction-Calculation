function design = loadBestMultiRadiusPhaseDesign(resultPath)
%LOADBESTMULTIRADIUSPHASEDESIGN 读取 DOE 最优多半径相位微岛设计。
%   design = loadBestMultiRadiusPhaseDesign() 默认从 ../ImageForShow 中读取
%   multi_radius_phase_doe_results.mat；若 MAT 不存在，则回退读取 CSV。
%   读取 results 表时，会先按工程约束筛选，再在合格候选中排序：
%       ZeroOrderRetention >= 0.95
%       RingVisibilityRatio <= 1.03
%       HazeRisk <= 0.12
%   若没有候选满足全部约束，则退回到 Objective 最小候选。
%
%   输出字段可直接用于 CustomApertureWithPhase_fraunhofer_fft.m：
%       design.phaseIslandRadiusList
%       design.phasePatternMode
%       design.phaseAspectRatio
%       design.phaseFillFactor
%       design.phaseMinDistanceFactor
%       design.phaseMinDistance
%       design.phaseSeed

    if nargin < 1 || strlength(string(resultPath)) == 0
        srcDir = fileparts(mfilename('fullpath'));
        resultDir = fullfile(srcDir, '..', 'ImageForShow');
        matPath = fullfile(resultDir, 'multi_radius_phase_doe_results.mat');
        csvPath = fullfile(resultDir, 'multi_radius_phase_doe_results.csv');
    else
        resultPath = char(resultPath);
        if isfolder(resultPath)
            matPath = fullfile(resultPath, 'multi_radius_phase_doe_results.mat');
            csvPath = fullfile(resultPath, 'multi_radius_phase_doe_results.csv');
        else
            [folderPath, ~, ext] = fileparts(resultPath);
            if strcmpi(ext, '.mat')
                matPath = resultPath;
                csvPath = fullfile(folderPath, 'multi_radius_phase_doe_results.csv');
            elseif strcmpi(ext, '.csv')
                csvPath = resultPath;
                matPath = fullfile(folderPath, 'multi_radius_phase_doe_results.mat');
            else
                error('不支持的 DOE 结果文件类型：%s', resultPath);
            end
        end
    end

    if exist(matPath, 'file')
        data = load(matPath, 'results', 'bestDesign', 'hazeRiskLimit', 'doePreset');
        if isfield(data, 'results') && ~isempty(data.results)
            [bestRow, selectionInfo] = selectBestDesignRow(data.results);
            design = rowToDesign(bestRow);
            design.selectionInfo = selectionInfo;
        elseif isfield(data, 'bestDesign')
            design = bestDesignToDesign(data.bestDesign);
            design.selectionInfo = "fallback_bestDesign_no_results_table";
        else
            error('MAT 文件中没有 results 或 bestDesign：%s', matPath);
        end
        design.sourcePath = matPath;
        if isfield(data, 'hazeRiskLimit')
            design.hazeRiskLimit = data.hazeRiskLimit;
        end
        if isfield(data, 'doePreset')
            design.doePreset = data.doePreset;
        end
        return
    end

    if exist(csvPath, 'file')
        results = readtable(csvPath);
        [bestRow, selectionInfo] = selectBestDesignRow(results);
        design = rowToDesign(bestRow);
        design.selectionInfo = selectionInfo;
        design.sourcePath = csvPath;
        return
    end

    error('未找到 DOE 结果文件：%s 或 %s', matPath, csvPath);
end

function [bestRow, selectionInfo] = selectBestDesignRow(results)
    %SELECTBESTDESIGNROW 先做约束筛选，再排序选择候选。
    zeroOrderRetentionMin = 0.95;
    ringVisibilityRatioMax = 1.03;
    hazeRiskMax = 0.12;

    candidateMask = true(height(results), 1);
    if any(strcmp(results.Properties.VariableNames, 'ZeroOrderRetention'))
        candidateMask = candidateMask & results.ZeroOrderRetention >= zeroOrderRetentionMin;
    end
    if any(strcmp(results.Properties.VariableNames, 'RingVisibilityRatio'))
        candidateMask = candidateMask & results.RingVisibilityRatio <= ringVisibilityRatioMax;
    end
    if any(strcmp(results.Properties.VariableNames, 'HazeRisk'))
        candidateMask = candidateMask & results.HazeRisk <= hazeRiskMax;
    end

    if any(candidateMask)
        candidates = results(candidateMask, :);
        candidates = sortCandidateRows(candidates);
        bestRow = candidates(1, :);
        selectionInfo = sprintf("constraint_filtered:%d/%d", height(candidates), height(results));
    else
        candidates = sortrows(results, 'Objective', 'ascend');
        bestRow = candidates(1, :);
        selectionInfo = sprintf("fallback_objective:no_constraint_candidate/%d", height(results));
    end
end

function rows = sortCandidateRows(rows)
    %SORTCANDIDATEROWS 合格候选优先按非零级削峰排序，再看综合目标。
    % 如果存在 EffectiveFillBin，则先按实际填充率分桶去重，避免 TargetFillFactor
    % 不同但实际生成图案相同或接近的候选重复支配排序。
    rows = keepBestPerEffectiveFillBin(rows);
    if any(strcmp(rows.Properties.VariableNames, 'NonZeroTopKPeakRatio'))
        if any(strcmp(rows.Properties.VariableNames, 'ActualFillFactor'))
            rows = sortrows(rows, {'NonZeroTopKPeakRatio', 'Objective', 'ActualFillFactor'}, ...
                {'ascend', 'ascend', 'ascend'});
        else
            rows = sortrows(rows, {'NonZeroTopKPeakRatio', 'Objective'}, {'ascend', 'ascend'});
        end
    else
        rows = sortrows(rows, 'Objective', 'ascend');
    end
end

function rows = keepBestPerEffectiveFillBin(rows)
    if ~any(strcmp(rows.Properties.VariableNames, 'EffectiveFillBin'))
        return
    end

    sortedRows = sortrows(rows, 'Objective', 'ascend');
    [~, firstIndex] = unique(sortedRows.EffectiveFillBin, 'stable');
    rows = sortedRows(sort(firstIndex), :);
end

function design = rowToDesign(row)
    radiusSet_um = parseRadiusSet(row.RadiusSet_um);
    design = struct();
    design.phasePatternMode = "circle";
    design.phaseAspectRatio = 1.0;
    if any(strcmp(row.Properties.VariableNames, 'PatternMode'))
        design.phasePatternMode = string(row.PatternMode);
    end
    if any(strcmp(row.Properties.VariableNames, 'AspectRatio'))
        design.phaseAspectRatio = row.AspectRatio;
    end
    design.radiusSet_um = radiusSet_um;
    design.phaseIslandRadiusList = radiusSet_um * 1e-6;
    if any(strcmp(row.Properties.VariableNames, 'PhaseHeight_nm'))
        design.phaseHeight = row.PhaseHeight_nm * 1e-9;
    end
    design.phaseFillFactor = row.TargetFillFactor;
    design.phaseMinDistanceFactor = row.MinDistanceFactor;
    design.phaseMinDistance = design.phaseMinDistanceFactor * 2 * ...
        computeDistanceRadius(design.phaseIslandRadiusList, design.phasePatternMode, design.phaseAspectRatio);
    design.phaseSeed = row.Seed;
    design.objective = row.Objective;

    optionalFields = {'DesignIndex', 'PatternMode', 'AspectRatio', 'PhaseHeight_nm', 'ActualFillFactor', ...
        'FillAchievementRatio', 'FillShortfall', 'EffectiveFillBin', ...
        'NonZeroTopKPeakRatio', 'TopKPeakRatio', 'PeakBackgroundRatioPenalty', 'ZeroOrderRetention', ...
        'NearRingEnergyRatio', 'NearRingEnergyPenalty', 'RingVisibility', ...
        'RingVisibilityRatio', 'RingVisibilityPenalty', 'HazeRisk', 'Objective'};
    for idx = 1:numel(optionalFields)
        fieldName = optionalFields{idx};
        if any(strcmp(row.Properties.VariableNames, fieldName))
            design.(fieldName) = row.(fieldName);
        end
    end
end

function design = bestDesignToDesign(bestDesign)
    design = struct();
    design.phasePatternMode = "circle";
    design.phaseAspectRatio = 1.0;
    if isfield(bestDesign, 'patternMode')
        design.phasePatternMode = bestDesign.patternMode;
    end
    if isfield(bestDesign, 'aspectRatio')
        design.phaseAspectRatio = bestDesign.aspectRatio;
    end
    design.radiusSet_um = bestDesign.radiusSet_um;
    design.phaseIslandRadiusList = bestDesign.radiusSet_um * 1e-6;
    if isfield(bestDesign, 'phaseHeight')
        design.phaseHeight = bestDesign.phaseHeight;
    end
    design.phaseFillFactor = bestDesign.fillFactor;
    design.phaseMinDistanceFactor = bestDesign.minDistanceFactor;
    design.phaseMinDistance = design.phaseMinDistanceFactor * 2 * ...
        computeDistanceRadius(design.phaseIslandRadiusList, design.phasePatternMode, design.phaseAspectRatio);
    design.phaseSeed = bestDesign.seed;
    design.ActualFillFactor = bestDesign.actualFillFactor;
    if isfield(bestDesign, 'metric')
        design.metric = bestDesign.metric;
        design.objective = bestDesign.metric.Objective;
    end
end

function distanceRadius = computeDistanceRadius(sizeList, patternMode, aspectRatio)
    switch string(patternMode)
        case {"random_ellipse", "blue_noise_ellipse"}
            distanceRadius = max(sizeList) / sqrt(aspectRatio);
        otherwise
            distanceRadius = max(sizeList);
    end
end

function radiusSet_um = parseRadiusSet(value)
    textValue = char(string(value));
    textValue = strrep(textValue, '[', ' ');
    textValue = strrep(textValue, ']', ' ');
    textValue = strrep(textValue, ',', ' ');
    radiusSet_um = sscanf(textValue, '%f').';
    if isempty(radiusSet_um)
        error('无法解析 RadiusSet_um：%s', char(string(value)));
    end
end
