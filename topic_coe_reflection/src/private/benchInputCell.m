function [cfg,native] = benchInputCell(cfg,stack,inputDir)
%benchInputCell 按原图元数据取得完整输入单元，不展开目标阵列。
    native=cfg.samplesPerPeriod;
    if isfield(stack.geometry,'inputSamplesPerPeriod') && stack.geometry.inputSamplesPerPeriod>0
        native=stack.geometry.inputSamplesPerPeriod;
    end
    first=stack.layers.(stack.layerNames(1)).file;
    files=dir(fullfile(inputDir,first+".*"));
    if isfile(fullfile(inputDir,first)), file=fullfile(inputDir,first);
    else
        assert(isscalar(files),'COEReflection:StackImage','基准首层请使用唯一图片文件名。');
        file=fullfile(files.folder,files.name);
    end
    info=imfinfo(file);
    periods=[info.Width,info.Height]/native;
    assert(all(mod(periods,1)==0) && all(mod(cfg.periods,periods)==0), ...
        'COEReflection:BenchmarkPeriodicity','目标须为原图完整Pitch单元的整数重复。');
    assert(stack.geometry.repeatToFit || isequal(cfg.periods,periods), ...
        'COEReflection:BenchmarkPeriodicity','周期基准要求显式repeatToFit或目标就是原图单元。');
    cfg.periods=periods;
end
