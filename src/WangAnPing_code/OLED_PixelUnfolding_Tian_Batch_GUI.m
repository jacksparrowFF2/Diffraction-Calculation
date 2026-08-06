function OLED_PixelUnfolding_Tian_Batch_GUI()
% 像素展开批量仿真 GUI —— 田字形双像素 + 视距MTF + 白平衡 + 评估
% 支持选择图片文件夹、像素模板A/B、输出文件夹
% 田字形排列：mod(i+j,2)==0 -> 模板A，否则 -> 模板B

%% 主窗口
fig = figure('Name', 'OLED像素展开仿真 (田字形双像素)', ...
    'NumberTitle', 'off', ...
    'Position', [100, 100, 950, 780], ...
    'MenuBar', 'none', ...
    'ToolBar', 'none', ...
    'Resize', 'off');
panel_bg = [1 1 1];

%% 左侧面板：图片选择
pnl_left = uipanel('Parent', fig, 'Title', '待仿真图片文件夹', ...
    'Position', [0.02, 0.45, 0.35, 0.53], 'BackgroundColor', panel_bg);
uicontrol('Parent', pnl_left, 'Style', 'pushbutton', ...
    'String', '选择图片文件夹', ...
    'Position', [10, 340, 120, 30], 'Callback', @cb_select_src_folder);
txt_src_folder = uicontrol('Parent', pnl_left, 'Style', 'text', ...
    'String', '未选择', 'Position', [140, 345, 180, 20], ...
    'BackgroundColor', panel_bg, 'HorizontalAlignment', 'left');
listbox_imgs = uicontrol('Parent', pnl_left, 'Style', 'listbox', ...
    'String', {}, 'Position', [10, 30, 310, 300], 'Max', 100, 'Min', 0);

%% 右侧面板：模板与参数
pnl_right = uipanel('Parent', fig, 'Title', '像素模板与仿真参数', ...
    'Position', [0.39, 0.45, 0.59, 0.53], 'BackgroundColor', panel_bg);
% 模板A
uicontrol('Parent', pnl_right, 'Style', 'text', ...
    'String', '模板 A:', 'Position', [10, 345, 60, 20], ...
    'BackgroundColor', panel_bg, 'HorizontalAlignment', 'left');
uicontrol('Parent', pnl_right, 'Style', 'pushbutton', ...
    'String', '选择文件', 'Position', [70, 345, 70, 25], ...
    'Callback', @cb_select_tmpl_A);
txt_tmpl_A = uicontrol('Parent', pnl_right, 'Style', 'text', ...
    'String', '未选择', 'Position', [150, 350, 250, 20], ...
    'BackgroundColor', panel_bg, 'HorizontalAlignment', 'left');
% 模板B
uicontrol('Parent', pnl_right, 'Style', 'text', ...
    'String', '模板 B:', 'Position', [10, 310, 60, 20], ...
    'BackgroundColor', panel_bg, 'HorizontalAlignment', 'left');
uicontrol('Parent', pnl_right, 'Style', 'pushbutton', ...
    'String', '选择文件', 'Position', [70, 310, 70, 25], ...
    'Callback', @cb_select_tmpl_B);
txt_tmpl_B = uicontrol('Parent', pnl_right, 'Style', 'text', ...
    'String', '未选择', 'Position', [150, 315, 250, 20], ...
    'BackgroundColor', panel_bg, 'HorizontalAlignment', 'left');
% 模板尺寸显示
txt_tmpl_size = uicontrol('Parent', pnl_right, 'Style', 'text', ...
    'String', '模板尺寸: 待加载', 'Position', [10, 280, 200, 20], ...
    'BackgroundColor', panel_bg, 'HorizontalAlignment', 'left');
% 像素间距与观看距离
uicontrol('Parent', pnl_right, 'Style', 'text', ...
    'String', '像素间距 (μm):', 'Position', [10, 250, 100, 20], ...
    'BackgroundColor', panel_bg, 'HorizontalAlignment', 'left');
edit_pitch = uicontrol('Parent', pnl_right, 'Style', 'edit', ...
    'String', '96', 'Position', [120, 250, 60, 25], 'BackgroundColor', 'white');
uicontrol('Parent', pnl_right, 'Style', 'text', ...
    'String', '观看距离 (m):', 'Position', [200, 250, 100, 20], ...
    'BackgroundColor', panel_bg, 'HorizontalAlignment', 'left');
edit_dist = uicontrol('Parent', pnl_right, 'Style', 'edit', ...
    'String', '0.45', 'Position', [310, 250, 60, 25], 'BackgroundColor', 'white');
% MTF模糊复选框
chk_mtf = uicontrol('Parent', pnl_right, 'Style', 'checkbox', ...
    'String', '启用 MTF 模糊', 'Position', [10, 220, 150, 25], ...
    'Value', 1, 'BackgroundColor', panel_bg);
% 最大边长
uicontrol('Parent', pnl_right, 'Style', 'text', ...
    'String', '最大边长 (px):', 'Position', [10, 190, 100, 20], ...
    'BackgroundColor', panel_bg, 'HorizontalAlignment', 'left');
edit_max_size = uicontrol('Parent', pnl_right, 'Style', 'edit', ...
    'String', '400', 'Position', [120, 190, 60, 25], 'BackgroundColor', 'white');
% 保存选项
chk_macro = uicontrol('Parent', pnl_right, 'Style', 'checkbox', ...
    'String', '保存微距视图', 'Position', [10, 160, 150, 25], ...
    'Value', 1, 'BackgroundColor', panel_bg);
chk_comp = uicontrol('Parent', pnl_right, 'Style', 'checkbox', ...
    'String', '保存对比图', 'Position', [170, 160, 150, 25], ...
    'Value', 1, 'BackgroundColor', panel_bg);

%% 底部面板：输出与评价
pnl_bottom = uipanel('Parent', fig, 'Title', '输出文件夹与评价', ...
    'Position', [0.02, 0.02, 0.96, 0.42], 'BackgroundColor', panel_bg);
uicontrol('Parent', pnl_bottom, 'Style', 'pushbutton', ...
    'String', '选择输出文件夹', 'Position', [10, 210, 120, 30], ...
    'Callback', @cb_select_output);
txt_out_folder = uicontrol('Parent', pnl_bottom, 'Style', 'text', ...
    'String', '未选择', 'Position', [140, 215, 300, 20], ...
    'BackgroundColor', panel_bg, 'HorizontalAlignment', 'left');
btn_run = uicontrol('Parent', pnl_bottom, 'Style', 'pushbutton', ...
    'String', '▶ 开始批量仿真', 'Position', [500, 200, 130, 45], ...
    'FontWeight', 'bold', 'FontSize', 12, ...
    'BackgroundColor', [0.3 0.7 0.3], 'ForegroundColor', 'white', ...
    'Callback', @cb_run);
progress_bar = uicontrol('Parent', pnl_bottom, 'Style', 'text', ...
    'String', '', 'Position', [10, 175, 850, 20], ...
    'BackgroundColor', [0.8 0.8 0.8]);
txt_status = uicontrol('Parent', pnl_bottom, 'Style', 'text', ...
    'String', '就绪', 'Position', [10, 150, 800, 20], ...
    'BackgroundColor', panel_bg, 'FontWeight', 'bold');
uicontrol('Parent', pnl_bottom, 'Style', 'text', ...
    'String', '评价结论 (平均值):', 'Position', [10, 125, 150, 20], ...
    'BackgroundColor', panel_bg, 'FontWeight', 'bold');
txt_eval = uicontrol('Parent', pnl_bottom, 'Style', 'text', ...
    'String', '等待仿真...', 'Position', [10, 5, 850, 115], ...
    'BackgroundColor', [1 1 1], 'HorizontalAlignment', 'left', ...
    'FontName', 'Consolas', 'FontSize', 9);

% 共享数据
handles = struct();
handles.src_folder = '';
handles.tmplA_path = '';
handles.tmplB_path = '';
handles.output_folder = '';
handles.listbox_imgs = listbox_imgs;
handles.txt_src_folder = txt_src_folder;
handles.txt_tmpl_A = txt_tmpl_A;
handles.txt_tmpl_B = txt_tmpl_B;
handles.txt_tmpl_size = txt_tmpl_size;
handles.edit_pitch = edit_pitch;
handles.edit_dist = edit_dist;
handles.chk_mtf = chk_mtf;
handles.edit_max_size = edit_max_size;
handles.chk_macro = chk_macro;
handles.chk_comp = chk_comp;
handles.txt_out_folder = txt_out_folder;
handles.btn_run = btn_run;
handles.progress_bar = progress_bar;
handles.txt_status = txt_status;
handles.txt_eval = txt_eval;
guidata(fig, handles);

%% 回调函数
    function cb_select_src_folder(~,~)
        handles = guidata(fig);
        folder = uigetdir(pwd, '选择包含图片的文件夹');
        if folder == 0, return; end
        handles.src_folder = folder;
        handles.txt_src_folder.String = folder;
        imgs = dir(fullfile(folder, '*.png'));
        if isempty(imgs)
            handles.listbox_imgs.String = {};
        else
            handles.listbox_imgs.String = {imgs.name};
        end
        guidata(fig, handles);
    end

    function cb_select_tmpl_A(~,~)
        handles = guidata(fig);
        [file, path] = uigetfile({'*.png;*.jpg;*.bmp', '图像文件'}, '选择模板 A');
        if isequal(file,0), return; end
        handles.tmplA_path = fullfile(path, file);
        handles.txt_tmpl_A.String = handles.tmplA_path;
        update_tmpl_size(handles);
        guidata(fig, handles);
    end

    function cb_select_tmpl_B(~,~)
        handles = guidata(fig);
        [file, path] = uigetfile({'*.png;*.jpg;*.bmp', '图像文件'}, '选择模板 B');
        if isequal(file,0), return; end
        handles.tmplB_path = fullfile(path, file);
        handles.txt_tmpl_B.String = handles.tmplB_path;
        update_tmpl_size(handles);
        guidata(fig, handles);
    end

    function update_tmpl_size(handles)
        if ~isempty(handles.tmplA_path) && ~isempty(handles.tmplB_path)
            try
                A = imread(handles.tmplA_path);
                B = imread(handles.tmplB_path);
                C = min(size(A,1), size(A,2));
                D = min(size(B,1), size(B,2));
                if C==D
                    handles.txt_tmpl_size.String = sprintf('模板尺寸: %d x %d', C, C);
                else
                    handles.txt_tmpl_size.String = '警告: A/B 尺寸不一致';
                end
            catch
                handles.txt_tmpl_size.String = '模板加载失败';
            end
        end
    end

    function cb_select_output(~,~)
        handles = guidata(fig);
        folder = uigetdir(pwd, '选择输出文件夹');
        if folder == 0, return; end
        handles.output_folder = folder;
        handles.txt_out_folder.String = folder;
        guidata(fig, handles);
    end

    function cb_run(~,~)
        handles = guidata(fig);
        if isempty(handles.src_folder)
            errordlg('请选择图片文件夹', '错误'); return;
        end
        img_list = handles.listbox_imgs.String;
        if isempty(img_list)
            errordlg('文件夹中没有 PNG 图片', '错误'); return;
        end
        if isempty(handles.tmplA_path) || isempty(handles.tmplB_path)
            errordlg('请选择模板 A 和 B', '错误'); return;
        end
        if isempty(handles.output_folder)
            errordlg('请选择输出文件夹', '错误'); return;
        end
        
        pitch_um = str2double(handles.edit_pitch.String);
        dist_m = str2double(handles.edit_dist.String);
        if isnan(pitch_um)||pitch_um<=0, errordlg('无效像素间距'); return; end
        if isnan(dist_m)||dist_m<=0, errordlg('无效观看距离'); return; end
        max_img_size = str2double(handles.edit_max_size.String);
        if isnan(max_img_size)||max_img_size<50, max_img_size=400; end
        enable_mtf = handles.chk_mtf.Value;
        
        % 加载模板
        try
            pxl_A = im2double(imread(handles.tmplA_path));
            pxl_B = im2double(imread(handles.tmplB_path));
        catch
            errordlg('模板加载失败', '错误'); return;
        end
        if size(pxl_A,3)==1, pxl_A=repmat(pxl_A,[1,1,3]); end
        if size(pxl_B,3)==1, pxl_B=repmat(pxl_B,[1,1,3]); end
        C = min(min(size(pxl_A,1),size(pxl_A,2)), min(size(pxl_B,1),size(pxl_B,2)));
        pxl_A = pxl_A(1:C, 1:C, :);
        pxl_B = pxl_B(1:C, 1:C, :);
        handles.txt_tmpl_size.String = sprintf('模板尺寸: %d x %d', C, C);
        
        % 计算平均有效透过率 (白平衡用)
        linA = srgb2linear(pxl_A);
        linB = srgb2linear(pxl_B);
        R_total = (mean(linA(:,:,1),'all') + mean(linB(:,:,1),'all'))/2;
        G_total = (mean(linA(:,:,2),'all') + mean(linB(:,:,2),'all'))/2;
        B_total = (mean(linA(:,:,3),'all') + mean(linB(:,:,3),'all'))/2;
        R_total = max(R_total,1e-6); G_total=max(G_total,1e-6); B_total=max(B_total,1e-6);
        gain = [1/R_total, 1/G_total, 1/B_total];
        
        set(handles.btn_run,'Enable','off');
        handles.txt_status.String = '处理中...';
        drawnow;
        
        n_imgs = length(img_list);
        results = []; valid=0;
        sum_ssim=0; sum_psnr=0; sum_bright=0; sum_contrast=0;
        sum_chroma=[0,0,0]; sum_jag=0; sum_sharp=0;
        
        for k = 1:n_imgs
            fname = img_list{k};
            handles.progress_bar.String = sprintf('处理 %d/%d: %s', k, n_imgs, fname);
            handles.progress_bar.BackgroundColor = [0.2 0.6 1.0];
            drawnow;
            
            img_path = fullfile(handles.src_folder, fname);
            fprintf('处理: %s ...', fname);
            try
                img_src = im2double(imread(img_path));
                if size(img_src,3)==1, img_src=repmat(img_src,[1,1,3]); end
                [H, W, ~] = size(img_src);
                if max(H,W) > max_img_size
                    scale = max_img_size / max(H,W);
                    img_src = imresize(img_src, scale);
                    [H, W, ~] = size(img_src);
                    fprintf(' (缩放至 %.0f%%)', scale*100);
                end
                
                % ---------- 像素展开 ----------
                img_lin = srgb2linear(img_src);
                % 信号扩展
                sig_R = kron(img_lin(:,:,1), ones(C,C));
                sig_G = kron(img_lin(:,:,2), ones(C,C));
                sig_B = kron(img_lin(:,:,3), ones(C,C));
                signal_exp = cat(3, sig_R, sig_G, sig_B);
                
                % 田字形模板平铺
                tile = zeros(H*C, W*C, 3);
                for i = 1:H
                    for j = 1:W
                        if mod(i+j,2)==0
                            blk = linA;
                        else
                            blk = linB;
                        end
                        r_idx = (i-1)*C+1;
                        c_idx = (j-1)*C+1;
                        tile(r_idx:r_idx+C-1, c_idx:c_idx+C-1, :) = blk;
                    end
                end
                
                % 调制
                highres_lin = signal_exp .* tile;
                
                % MTF模糊 (可选)
                if enable_mtf
                    sigma_pix = (dist_m * tan(1/60*pi/180)) / (pitch_um*1e-6);
                    sigma_pix = max(sigma_pix, 0.5);
                    highres_lin = imgaussfilt(highres_lin, sigma_pix);
                end
                
                % 视觉积分
                lowres_lin = block_average(highres_lin, C);
                
                % 白平衡后补偿
                lowres_comp = lowres_lin;
                for c = 1:3
                    lowres_comp(:,:,c) = lowres_lin(:,:,c) * gain(c);
                end
                lowres_comp = max(0, min(1, lowres_comp));
                
                % 转 sRGB
                sim_img = linear2srgb(lowres_comp);
                highres_srgb = linear2srgb(highres_lin);
                
                % 保存
                [~,name,~] = fileparts(fname);
                out_sub = fullfile(handles.output_folder, name);
                if ~exist(out_sub,'dir'), mkdir(out_sub); end
                imwrite(img_src, fullfile(out_sub,'original.png'));
                imwrite(sim_img, fullfile(out_sub,'simulated.png'));
                if handles.chk_macro.Value
                    imwrite(highres_srgb, fullfile(out_sub,'macro_view.png'));
                end
                if handles.chk_comp.Value
                    comp = [img_src, sim_img];
                    imwrite(comp, fullfile(out_sub,'comparison.png'));
                end
                
                % 评估
                orig_gray = rgb2gray(img_src);
                sim_gray  = rgb2gray(sim_img);
                ssim_val = ssim(sim_gray, orig_gray);
                psnr_val = psnr(sim_img, img_src);
                lum_orig = mean(orig_gray(:));
                lum_sim  = mean(sim_gray(:));
                bright_ret = (lum_sim / lum_orig) * 100;
                orig_contrast = std2(orig_gray);
                sim_contrast = std2(sim_gray);
                contrast_ret = (sim_contrast / max(orig_contrast,1e-6)) * 100;
                chroma_err = squeeze(mean(abs(sim_img - img_src), [1,2]))';
                if numel(chroma_err)<3, chroma_err(end+1:3)=0; end
                grad_orig = imgradient(orig_gray);
                grad_sim  = imgradient(sim_gray);
                jag_ratio = mean(grad_sim(:)) / mean(grad_orig(:));
                if jag_ratio<1.1, jagginess='低'; elseif jag_ratio<1.4, jagginess='中'; else, jagginess='高'; end
                lap_orig = mean(abs(imfilter(orig_gray, fspecial('laplacian'))), 'all');
                lap_sim  = mean(abs(imfilter(sim_gray, fspecial('laplacian'))), 'all');
                sharp_ratio = lap_sim / lap_orig;
                if sharp_ratio>0.9, sharpness='极高'; elseif sharp_ratio>0.7, sharpness='高';
                elseif sharp_ratio>0.5, sharpness='中'; else, sharpness='低'; end
                if ssim_val>0.95, fidelity='极高'; elseif ssim_val>0.9, fidelity='高';
                elseif ssim_val>0.7, fidelity='中'; else, fidelity='低'; end
                
                valid = valid+1;
                results(valid).filename = fname;
                results(valid).width = W; results(valid).height = H;
                results(valid).ssim = ssim_val;
                results(valid).psnr = psnr_val;
                results(valid).bright_ret = bright_ret;
                results(valid).contrast_ret = contrast_ret;
                results(valid).jag_ratio = jag_ratio;
                results(valid).sharp_ratio = sharp_ratio;
                results(valid).chroma_R = chroma_err(1);
                results(valid).chroma_G = chroma_err(2);
                results(valid).chroma_B = chroma_err(3);
                results(valid).fidelity = fidelity;
                results(valid).jagginess = jagginess;
                results(valid).sharpness = sharpness;
                
                sum_ssim = sum_ssim+ssim_val;
                sum_psnr = sum_psnr+psnr_val;
                sum_bright = sum_bright+bright_ret;
                sum_contrast = sum_contrast+contrast_ret;
                sum_chroma = sum_chroma+chroma_err(1:3);
                sum_jag = sum_jag+jag_ratio;
                sum_sharp = sum_sharp+sharp_ratio;
                
                fprintf(' 成功 (SSIM: %.3f)\n', ssim_val);
            catch ME
                fprintf(' 失败: %s\n', ME.message);
            end
        end
        
        if valid>0
            avg_ssim=sum_ssim/valid; avg_psnr=sum_psnr/valid;
            avg_bright=sum_bright/valid; avg_contrast=sum_contrast/valid;
            avg_chroma=sum_chroma/valid; avg_jag=sum_jag/valid; avg_sharp=sum_sharp/valid;
            fids = {results.fidelity};
            high=sum(strcmp(fids,'极高')); good=sum(strcmp(fids,'高'));
            mid=sum(strcmp(fids,'中')); low=sum(strcmp(fids,'低'));
            eval_str = sprintf([ ...
                '平均 SSIM: %.4f | 平均 PSNR: %.2f dB\n', ...
                '亮度保持率: %.2f%% | 对比度保持率: %.2f%%\n', ...
                '色度误差 (R,G,B): %.4f, %.4f, %.4f\n', ...
                '锯齿度比率: %.3f | 锐利度比率: %.3f\n', ...
                '保真度分布: 极高 %d, 高 %d, 中 %d, 低 %d'], ...
                avg_ssim, avg_psnr, avg_bright, avg_contrast, ...
                avg_chroma(1), avg_chroma(2), avg_chroma(3), ...
                avg_jag, avg_sharp, high, good, mid, low);
            handles.txt_eval.String = eval_str;
            fprintf('\n=== 评价结论 ===\n%s\n', eval_str);
            export_excel(results, handles.output_folder, struct('C',C,'R_total',R_total,...
                'G_total',G_total,'B_total',B_total,'pitch_um',pitch_um,'dist_m',dist_m,'mtf',enable_mtf));
            handles.txt_status.String = sprintf('完成 %d/%d 张，报告已生成', valid, n_imgs);
        else
            handles.txt_eval.String = '无成功处理图片';
            handles.txt_status.String = '仿真失败';
        end
        handles.progress_bar.String = '完成';
        handles.progress_bar.BackgroundColor = [0.3 0.8 0.3];
        set(handles.btn_run,'Enable','on');
        guidata(fig, handles);
        drawnow;
    end
end

%% 辅助函数
function lin = srgb2linear(srgb)
    lin = zeros(size(srgb));
    mask = srgb <= 0.04045;
    lin(mask) = srgb(mask) / 12.92;
    lin(~mask) = ((srgb(~mask) + 0.055) / 1.055) .^ 2.4;
end

function srgb = linear2srgb(lin)
    srgb = zeros(size(lin));
    mask = lin <= 0.0031308;
    srgb(mask) = lin(mask) * 12.92;
    srgb(~mask) = 1.055 * (lin(~mask) .^ (1/2.4)) - 0.055;
end

function J = block_average(I, factor)
    [H, W, ~] = size(I);
    H_new = floor(H / factor);
    W_new = floor(W / factor);
    I = I(1:H_new*factor, 1:W_new*factor, :);
    J = zeros(H_new, W_new, 3);
    for c = 1:3
        temp = reshape(I(:,:,c), factor, H_new, factor, W_new);
        J(:,:,c) = squeeze(mean(mean(temp, 1), 3));
    end
end

function export_excel(results, out_folder, info)
    if isempty(results), return; end
    n = numel(results);
    filename = cell(n,1); width=zeros(n,1); height=zeros(n,1);
    ssim_val=zeros(n,1); psnr_val=zeros(n,1); bright=zeros(n,1); contrast=zeros(n,1);
    chroma=zeros(n,3); jag=zeros(n,1); sharp=zeros(n,1);
    fidelity=cell(n,1); jagginess=cell(n,1); sharpness=cell(n,1);
    for k=1:n
        r = results(k);
        filename{k}=r.filename; width(k)=r.width; height(k)=r.height;
        ssim_val(k)=r.ssim; psnr_val(k)=r.psnr;
        bright(k)=r.bright_ret; contrast(k)=r.contrast_ret;
        chroma(k,:)=[r.chroma_R, r.chroma_G, r.chroma_B];
        jag(k)=r.jag_ratio; sharp(k)=r.sharp_ratio;
        fidelity{k}=r.fidelity; jagginess{k}=r.jagginess; sharpness{k}=r.sharpness;
    end
    T = table(filename, width, height, ssim_val, psnr_val, bright, contrast, ...
        chroma(:,1), chroma(:,2), chroma(:,3), jag, sharp, fidelity, jagginess, sharpness, ...
        'VariableNames', {'文件名','宽度','高度','SSIM','PSNR_dB','亮度保持率_%','对比度保持率_%', ...
        '色度误差_R','色度误差_G','色度误差_B','锯齿度比率','锐利度比率', ...
        '保真度','锯齿度定性','锐利度定性'});
    T_param = cell2table({
        '模板尺寸', info.C;
        'R透过率', info.R_total;
        'G透过率', info.G_total;
        'B透过率', info.B_total;
        '像素间距(μm)', info.pitch_um;
        '观看距离(m)', info.dist_m;
        'MTF模糊', info.mtf;
        }, 'VariableNames', {'参数','值'});
    fullname = fullfile(out_folder, '仿真报告.xlsx');
    writetable(T, fullname, 'Sheet', '评估结果');
    writetable(T_param, fullname, 'Sheet', '仿真参数');
    fprintf('Excel 报告已保存至: %s\n', fullname);
end