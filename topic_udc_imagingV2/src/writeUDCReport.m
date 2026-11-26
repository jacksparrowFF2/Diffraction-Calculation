function reportPath = writeUDCReport(projectOutputDir,modes)
%writeUDCReport 汇总已完成模式的 CSV、配置和图像，生成可离线打开的 HTML。
    arguments
        projectOutputDir (1,1) string
        modes (1,:) string {mustBeNonempty}
    end
    modes = coreCanonicalUDCModeNames(modes);
    [~,projectName] = fileparts(projectOutputDir);
    heading = escapeHtml(projectName)+" · UDC 成像报告";
    sections = strings(numel(modes),1);
    navigation = sections;
    for k = 1:numel(modes)
        mode = modes(k);
        assert(isvarname(char(mode)), 'UDCImagingV2:ReportMode','无效的透射模式名：%s',mode);
        label = mode;
        if mode=="specifiedTransmittance", label = "给定功率透过率"; end
        if mode=="nk", label = "薄膜 n、k 与厚度"; end
        label = escapeHtml(label);
        navigation(k) = "<a href=""#"+mode+""">"+label+"</a>";
        modeOutputDir = coreFindUDCModeResults(projectOutputDir,mode);
        sections(k) = modeSection(modeOutputDir,mode,label);
    end

    style = [
        "html{background:#f4f6fa;color:#17263c;font:16px/1.6 system-ui,sans-serif}"
        "body{max-width:1120px;margin:32px auto;padding:0 20px}h1,h2,h3{line-height:1.2}"
        "nav{display:flex;gap:18px;flex-wrap:wrap}a{color:#155e9c}"
        "section{margin:26px 0;padding:24px;background:white;border:1px solid #dbe3ec;border-radius:12px}"
        ".scroll{overflow-x:auto}table{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums}"
        "caption{text-align:left;font-weight:650;margin:14px 0 8px}"
        "th,td{padding:7px 10px;text-align:right;white-space:nowrap;border-bottom:1px solid #dbe3ec}"
        "th{background:#eef3f9}th:first-child,td:first-child{text-align:left}"
        ".gallery{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:14px}"
        ".comparison{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px}"
        "figure{margin:15px 0}img{width:100%;height:auto;display:block;border-radius:5px}"
        "figcaption{margin-top:5px;color:#40566e}"
        "@media print{html{background:white}section{padding:0;border:none}figure,tr{break-inside:avoid}}"
    ];
    html = "<!doctype html><html lang=""zh-CN""><head><meta charset=""utf-8"">"+ ...
        "<meta name=""viewport"" content=""width=device-width, initial-scale=1"">"+ ...
        "<title>"+heading+"</title><style>"+join(style,newline)+"</style></head>"+ ...
        "<body><h1>"+heading+"</h1><nav aria-label=""透射模式"">"+join(navigation,"")+ ...
        "</nav><main>"+join(sections,"")+"</main></body></html>";

    reportPath = fullfile(projectOutputDir,'udc_report.html');
    temporaryPath = string(tempname(projectOutputDir));
    cleanTemporary = onCleanup(@() removeTemporary(temporaryPath));
    file = openFile(temporaryPath,'w');
    closeFile = onCleanup(@() fclose(file));
    fprintf(file,'%s',html);
    [message,errorNumber] = ferror(file);
    assert(errorNumber==0,'UDCImagingV2:ReportWrite','写入报告失败：%s',message);
    clear closeFile
    [ok,message] = movefile(temporaryPath,reportPath,'f');
    assert(ok,'UDCImagingV2:ReportWrite','保存报告失败：%s',message);
end

function html = modeSection(folder,mode,label)
    data = load(fullfile(folder,'udc_results.mat'),'cfg');
    cfg = data.cfg;
    comparisonImage = "udc_imaging_comparison.png";
    if ~isfile(fullfile(folder,comparisonImage))
        existing = load(fullfile(folder,'udc_results.mat'),'optics');
        plotImagingComparison(existing.optics,folder);
        clear existing
    end
    spanMm = (cfg.outputSize-1)*cfg.sensorPitch*1e3;
    setup = sprintf(['阵列 %d×%d Pitch，入瞳直径 %.3g mm，焦距 %.3g mm，像元 %.3g μm。' ...
        'PSF 网格 %d×%d，像元中心跨度 %.3g mm。波长 %.0f–%.0f nm，共 %d 点。'], ...
        cfg.periods(1),cfg.periods(2),2e3*cfg.pupilRadius,1e3*cfg.focalLength, ...
        1e6*cfg.sensorPitch,cfg.outputSize,cfg.outputSize,spanMm, ...
        min(cfg.wavelengths)*1e9,max(cfg.wavelengths)*1e9,numel(cfg.wavelengths));
    if isfield(cfg,'reusedFromProject')
        setup = string(setup) + sprintf(' 复用 %s 的衍射解。理想中性层强度透过率为 %.3g%%，所有波长相同。', ...
            cfg.reusedFromProject,100*cfg.neutralTransmittance);
    end
    channel = readtable(fullfile(folder,'udc_channel_metrics.csv'),'TextType','string');
    mtf = readtable(fullfile(folder,'udc_mtf_directional_summary.csv'),'TextType','string');
    tail = readtable(fullfile(folder,'udc_psf_tail.csv'),'TextType','string');
    mtfBand = numberText(mtf.BandStartCyclesPerMm)+"–"+numberText(mtf.BandEndCyclesPerMm);
    directions = replace(mtf.Direction,["Horizontal","Vertical","Diagonal45","Diagonal135"], ...
        ["水平","垂直","对角 45°","对角 135°"]);

    html = "<section id="""+mode+"""><h2>"+label+"</h2><p>"+escapeHtml(setup)+"</p>"+ ...
        "<h3>同曝光成像对比</h3><div class=""comparison"">"+ ...
        figureHtml(folder,"udc_clear.png","无遮挡 · 仅入瞳与传感器像元积分")+ ...
        figureHtml(folder,"udc_fixed_exposure.png","含 UDC Mask · 同一曝光")+"</div>"+ ...
        figureHtml(folder,comparisonImage,"无遮挡与含 Mask 的 PSF 全范围、中心放大和 MTF")+ ...
        htmlTable("进光量与偏色", ...
        ["通道","总透过率","窗口捕获透过率","相对无遮挡","窗口内／总透射能量"], ...
        [channel.Channel,percentText(channel.Throughput),percentText(channel.CapturedThroughput), ...
        percentText(channel.RelativeToClear),percentText(channel.FractionOfTransmittedCaptured)])+ ...
        htmlTable("细节与方向性", ...
        ["通道","方向","MTF50 (cycles/mm)","无遮挡 MTF50","评价频段 (cycles/mm)","最低 MTF","无遮挡最低 MTF"], ...
        [mtf.Channel,directions,numberText(mtf.MTF50CyclesPerMm),numberText(mtf.ClearMTF50CyclesPerMm), ...
        mtfBand,numberText(mtf.BandMinMTF),numberText(mtf.ClearBandMinMTF)])+ ...
        htmlTable("光晕与拖尾", ...
        ["通道","半径 (μm)","半径外能量占比","无遮挡占比","尾部绝对透过率","无遮挡尾部透过率"], ...
        [tail.Channel,numberText(tail.RadiusUm),percentText(tail.TailFractionWithinWindow), ...
        percentText(tail.ClearTailFractionWithinWindow),percentText(tail.TailThroughput), ...
        percentText(tail.ClearTailThroughput)])+ ...
        "<p>尾部指标限于计算窗口。横线表示当前网格无法取得有效指标。</p>"+ ...
        figureHtml(folder,"udc_diffraction_diagnostics.png","有效入瞳、RGB 衍射强度、对数 PSF 与二维 MTF")+ ...
        figureHtml(folder,"udc_hdr_point.png","同曝光 HDR 点光源对照")+ ...
        "<div class=""gallery"">"+figureHtml(folder,"udc_brightness_matched.png","亮度匹配成像")+"</div></section>";
end

function text = numberText(values)
    text = compose("%.4g",values);
    text(~isfinite(values)) = "—";
end

function text = percentText(values)
    text = compose("%.2f%%",100*values);
    text(~isfinite(values)) = "—";
end

function html = htmlTable(label,headers,values)
    cells = "<td>"+escapeHtml(values)+"</td>";
    body = "<tr>"+join(cells,"",2)+"</tr>";
    html = "<div class=""scroll""><table><caption>"+label+"</caption><thead><tr>"+ ...
        join("<th scope=""col"">"+escapeHtml(headers)+"</th>","")+ ...
        "</tr></thead><tbody>"+join(body,"")+"</tbody></table></div>";
end

function html = figureHtml(folder,name,label)
    file = openFile(fullfile(folder,name),'r');
    closeFile = onCleanup(@() fclose(file));
    encoded = string(matlab.net.base64encode(fread(file,Inf,'*uint8')));
    label = escapeHtml(label);
    html = "<figure><img alt="""+label+""" src=""data:image/png;base64,"+encoded+ ...
        """><figcaption>"+label+"</figcaption></figure>";
end

function file = openFile(path,mode)
    [file,message] = fopen(path,mode,'n','UTF-8');
    assert(file>=0,'UDCImagingV2:ReportFile','不能打开 %s：%s',path,message);
end

function text = escapeHtml(text)
    text = replace(string(text),"&","&amp;");
    text = replace(text,["<",">",string(char(34))],["&lt;","&gt;","&quot;"]);
end

function removeTemporary(path)
    if isfile(path), delete(path); end
end
