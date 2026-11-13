function addCoeAlgorithmFixture(testCase, root)
%addCoeAlgorithmFixture 在临时副本中保留原有数值算法回归测试。
% 生产模块的 private 不加入路径；真实目录解析由 ModuleIsolationTest 验证。
    temporary = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
    accessDir = fullfile(temporary.Folder, 'src', 'algorithmAccess');
    mkdir(accessDir);
    module = fullfile(root, 'topic_coe_reflection');
    copyfile(fullfile(module, 'src', 'private', '*.m'), accessDir);
    % 保持算法文件到 Input 的相对层级，供 CIE 默认数据定位使用。
    copyfile(fullfile(module, 'Input'), fullfile(temporary.Folder, 'Input'));
    testCase.applyFixture(matlab.unittest.fixtures.PathFixture(accessDir));
end
