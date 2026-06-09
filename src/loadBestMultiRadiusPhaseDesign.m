function design = loadBestMultiRadiusPhaseDesign(resultPath)
%LOADBESTMULTIRADIUSPHASEDESIGN 读取 DOE 最优多半径相位微岛设计。
%   design = loadBestMultiRadiusPhaseDesign() 默认从 ../ImageForShow 中读取
%   multi_radius_phase_doe_results.mat；若 MAT 不存在，则回退读取 CSV。
%
%   输出字段可直接用于 CustomApertureWithPhase_fraunhofer_fft.m：
%       design.phaseIslandRadiusList
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
            results = sortrows(data.results, 'Objective', 'ascend');
            bestRow = results(1, :);
            design = rowToDesign(bestRow);
        elseif isfield(data, 'bestDesign')
            design = bestDesignToDesign(data.bestDesign);
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
        results = sortrows(results, 'Objective', 'ascend');
        design = rowToDesign(results(1, :));
        design.sourcePath = csvPath;
        return
    end

    error('未找到 DOE 结果文件：%s 或 %s', matPath, csvPath);
end

function design = rowToDesign(row)
    radiusSet_um = parseRadiusSet(row.RadiusSet_um);
    design = struct();
    design.radiusSet_um = radiusSet_um;
    design.phaseIslandRadiusList = radiusSet_um * 1e-6;
    if any(strcmp(row.Properties.VariableNames, 'PhaseHeight_nm'))
        design.phaseHeight = row.PhaseHeight_nm * 1e-9;
    end
    design.phaseFillFactor = row.TargetFillFactor;
    design.phaseMinDistanceFactor = row.MinDistanceFactor;
    design.phaseMinDistance = design.phaseMinDistanceFactor * 2 * max(design.phaseIslandRadiusList);
    design.phaseSeed = row.Seed;
    design.objective = row.Objective;

    optionalFields = {'DesignIndex', 'PhaseHeight_nm', 'ActualFillFactor', 'NonZeroTopKPeakRatio', ...
        'TopKPeakRatio', 'PeakBackgroundRatioPenalty', 'ZeroOrderRetention', ...
        'NearRingEnergyRatio', 'RingVisibility', 'RingVisibilityRatio', 'HazeRisk'};
    for idx = 1:numel(optionalFields)
        fieldName = optionalFields{idx};
        if any(strcmp(row.Properties.VariableNames, fieldName))
            design.(fieldName) = row.(fieldName);
        end
    end
end

function design = bestDesignToDesign(bestDesign)
    design = struct();
    design.radiusSet_um = bestDesign.radiusSet_um;
    design.phaseIslandRadiusList = bestDesign.radiusSet_um * 1e-6;
    if isfield(bestDesign, 'phaseHeight')
        design.phaseHeight = bestDesign.phaseHeight;
    end
    design.phaseFillFactor = bestDesign.fillFactor;
    design.phaseMinDistanceFactor = bestDesign.minDistanceFactor;
    design.phaseMinDistance = design.phaseMinDistanceFactor * 2 * max(design.phaseIslandRadiusList);
    design.phaseSeed = bestDesign.seed;
    design.ActualFillFactor = bestDesign.actualFillFactor;
    if isfield(bestDesign, 'metric')
        design.metric = bestDesign.metric;
        design.objective = bestDesign.metric.Objective;
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
