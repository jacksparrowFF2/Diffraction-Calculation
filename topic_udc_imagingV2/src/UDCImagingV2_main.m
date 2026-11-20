% 选择项目和模式。留空时运行项目配置中的全部模式。
projectName = "case01";
transmissionModes = strings(1,0);
[optics,images,reportPath] = runUDCProject(projectName,transmissionModes);
