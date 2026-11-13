function createCoeImageExamples(cfg, stack, inputRoot)
%createCoeImageExamples 导出三套同几何、不同图片格式的完整 Diamond 输入。
% 仅用于开发阶段预生成示例；用户运行主脚本时不调用，也不提供像素调整参数。
    stack.geometry.source = "diamond";
    [~, geometry] = prepareCoeStack(cfg, stack);
    fields = ["anodeMask", "bmMask", "cfRMask", "cfGMask", "cfBMask"];
    stems = ["anode", "BM", "CF_R", "CF_G", "CF_B"];
    extensions = ["png", "bmp", "jpg"];
    for j = 1:numel(extensions)
        project = "Case" + j;
        folder = fullfile(inputRoot, project);
        if ~isfolder(folder)
            mkdir(folder);
        end
        for k = 1:numel(fields)
            image = uint8(geometry.(fields(k)))*255;
            file = fullfile(folder, stems(k) + "." + extensions(j));
            if extensions(j) == "jpg"
                imwrite(image, file, 'Quality', 100);
            else
                imwrite(image, file);
            end
        end
        file = fopen(fullfile(folder, 'README.md'), 'w', 'n', 'UTF-8');
        assert(file >= 0, 'COEReflection:ExampleWrite', '无法写入示例说明。');
        cleanup = onCleanup(@() fclose(file));
        fprintf(file, '# %s：Diamond %s 示例图片集\n\n', project, upper(extensions(j)));
        fprintf(file, '这是演示几何，不代表实际产品。三套项目几何相同，仅编码格式不同。\n\n');
        fprintf(file, '- 主脚本设置 `projectName="%s"`、`stack.geometry.source="images"`。\n', project);
        fprintf(file, '- Pitch：%.6g μm；periods：[%d,%d]；samplesPerPeriod：%d。\n', ...
            cfg.pitch*1e6, cfg.periods, cfg.samplesPerPeriod);
        fprintf(file, '- 图片尺寸：%d 行 × %d 列；物面：%.6g × %.6g μm。\n', ...
            size(geometry.anodeMask), cfg.periods*cfg.pitch*1e6);
        fprintf(file, '- anode 白色=阳极；BM 白色=实心 BM（bmImageIsOpening=false）；CF_R/G/B 白色=对应膜。\n');
        fprintf(file, '- 阈值 0.5；图层文件名可不写扩展名，自动匹配本项目格式。\n');
        fprintf(file, '- 各项目 config.txt 自动提供采样、材料、观察角域和显示参数；不包含像素形状参数。\n');
        fprintf(file, '- 参考图实线 PDL 作为阳极，虚线 Mask 不作为 BM；R/B 交错、G 在四邻中心且左右交替倾斜，R:G:B 数量比 1:2:1。\n');
        fprintf(file, '- JPG 为有损格式，推荐实际二值版图优先使用 PNG/BMP；此示例使用 Quality=100。\n');
        fprintf(file, '- 材料仍取主脚本示例常数；若提供 nk CSV，请置于本项目 Materials/ 下并填写相对路径。\n');
        clear cleanup
    end
end
