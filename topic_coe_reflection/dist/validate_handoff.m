%% 验证解压后的 COE 迁移包（不依赖 Codex/MCP）
clearvars;
clc;
close all;
bundleRoot = fileparts(mfilename('fullpath'));
moduleSource = fullfile(bundleRoot, 'topic_coe_reflection', 'src');
assert(isfolder(moduleSource), 'COEReflection:BundleLayout', ...
    '请在完整迁移包根目录运行 validate_handoff.m。');
addpath(moduleSource);
assert(strcmpi(fileparts(which('runCOEReflection')), moduleSource), ...
    'COEReflection:BundlePath', 'MATLAB 未解析到当前迁移包的源码，请检查路径。');
suite = matlab.unittest.TestSuite.fromFolder(fullfile(bundleRoot, 'tests'));
assert(~isempty(suite), 'COEReflection:BundleTests', '未找到迁移包测试。');
validation = run(suite);
disp(table(validation));
assert(all([validation.Passed]), 'COEReflection:Validation', '迁移验证未全部通过，请查看失败信息。');
fprintf('迁移包验证通过。下一步运行 topic_coe_reflection/src/COEReflection_main.m。\n');
