%% F192 教学练习：只读 CSV 与小数组，不运行正式传播，不改写输入或输出
% 在编辑器打开此文件，逐节“运行节”。变量统一 lab 前缀，便于最后清理。
% 配套说明：F192_MATLAB教学手册.html。MATLAB R2024a 已核对。
labDocs = fileparts(mfilename('fullpath')); % 当前脚本所在 docs
labModule = fileparts(labDocs);             % 上一级模块目录
labRoot = fullfile(labModule,'Output','MaterialSpectra_20261117_nkEnabled');
assert(isfolder(labRoot),'请在完整项目中运行本练习，结果目录不存在。');

%% 1. MATLAB 数组基础（以下全为教学例子，不是正式光学结果）
labA = [1 2 3; 4 5 6]; % 两行三列；空格分列，分号分行
assert(isequal(size(labA),[2 3]));
disp(labA(:,2));        % 第2列：[2;5]
disp(labA(:));          % 按列展开：[1;4;2;5;3;6]
labCube = reshape(1:24,2,3,4); % 四个波长页，每页2×3
labPixel = squeeze(labCube(2,3,:)); % 4×1：[6;12;18;24]
assert(isequal(labPixel,[6;12;18;24]));
labWeights = [0.1;0.2;0.3;0.4]; % 4×1，权重和为1
labIntegrated = reshape(reshape(labCube,[],4)*labWeights,2,3);
assert(abs(labIntegrated(2,3)-18)<1e-12);
labR = sqrt(0.25)*exp(1i*pi); % 复振幅；模平方才是强度反射率
assert(abs(abs(labR)^2-0.25)<1e-12);

%% 2. 正式结果：由逐波长 CSV 与独立权重重算宽带指标
labMetrics = readtable(fullfile(labRoot,'F192','spectral_metrics.csv'));
labW = readtable(fullfile(labRoot,'spectrum_weights.csv'));
assert(isequal(labMetrics.WavelengthNm,labW.WavelengthNm));
assert(abs(sum(labW.EnergyWeight)-1)<1e-12);
assert(abs(sum(labW.PhotopicWeight)-1)<1e-12);
labTotal = dot(labMetrics.Throughput,labW.EnergyWeight);
labPhoto = dot(labMetrics.Throughput,labW.PhotopicWeight);
labWindow = dot(labMetrics.CapturedThroughput,labW.EnergyWeight);
% ZeroOrderFraction 的分母是该波长总反射，须先乘回 Throughput
labZero = dot(labMetrics.Throughput.*labMetrics.ZeroOrderFraction,labW.EnergyWeight);
labOutside = dot(labMetrics.Throughput.*labMetrics.OutsideAngleFraction,labW.EnergyWeight);
labCone = labTotal-labOutside;
labSummary = readtable(fullfile(labRoot,'comparison_summary.csv'));
labRow = labSummary(string(labSummary.Project)=="F192",:);
assert(abs(100*labTotal-labRow.RadiantTotalPercent)<1e-10);
assert(abs(100*labPhoto-labRow.PhotopicTotalPercent)<1e-10);
assert(abs(100*labWindow-labRow.RadiantWindowPercent)<1e-10);
assert(abs(100*labZero-labRow.RadiantZeroPercent)<1e-10);
assert(abs(100*labCone-labRow.RadiantCone2p5Percent)<1e-10);
fprintf('总反射 %.9f%%；窗口 %.9f%%；零级 %.9f%%；2.5度内 %.9f%%\n', ...
    100*labTotal,100*labWindow,100*labZero,100*labCone);
fprintf('窗口捕获/总反射 %.9f%%\n',100*labWindow/labTotal);

%% 3. 重建 XYZ 权重（41×3）；这里只处理 CSV，不加载能量立方体
labCMF = [labW.CMF_XYZ_1,labW.CMF_XYZ_2,labW.CMF_XYZ_3];
labRaw = labW.D65SPD.*labW.QuadratureNm; % SPD×梯形权重
labXYZWeights = labRaw.*labCMF/sum(labRaw.*labCMF(:,2));
assert(isequal(size(labXYZWeights),[41 3]));
assert(max(abs(labXYZWeights(:,2)-labW.PhotopicWeight))<1e-12);

%% 4. 单个假设像素的完整映射（教学数值，不是 F192 像素）
labXYZd = [120 100 80];
labGain = 0.0006609544398142476;
labMatrix = [3.2404542 -1.5371385 -0.4985314; ...
    -0.9692660 1.8760108 0.0415560;0.0556434 -0.2040259 1.0572252];
labYe = labXYZd(2)*labGain;
labYt = (log1p(60*labYe)/log(61))^0.72;
labScale = 1;
if labYe>0, labScale=labYt/labYe; end
labLinear = (labXYZd*labMatrix.')*labGain;
labEnhanced = labLinear*labScale;
labEnhanced = labEnhanced/max(1,max(labEnhanced));
labClipped = min(max(labEnhanced,0),1);
labEncoded = 12.92*labClipped;
labNonlinear = labClipped>0.0031308;
labEncoded(labNonlinear)=1.055*labClipped(labNonlinear).^(1/2.4)-0.055;
disp(table(labYe,labYt,labScale)); disp(labEncoded);

%% 5. 物理尺度：单位始终先用 m 和 rad
labPitch = 59.06e-6;
labDx = labPitch/64;
labLength = 4*labPitch;
labArea = labLength^2;
labTheta = 550e-9/labPitch;
labMillimetres = 0.25*labTheta*1e3;
fprintf('dx %.8f um；物面 %.5f um；名义间隔 %.8f deg / %.8f mm\n', ...
    labDx*1e6,labLength*1e6,rad2deg(labTheta),labMillimetres);
assert(abs(labArea-5.58093376e-8)<1e-20);

%% 6. 可选看图：取消下面注释，使用基础 image，无需 imshow
% labPNG = imread(fullfile(labRoot,'F192','rgb_angle_enhanced_shared.png'));
% figure; image(labPNG); axis image off;

%% 清理（手动选择运行；不清除用户的其它变量）
% clear lab*
