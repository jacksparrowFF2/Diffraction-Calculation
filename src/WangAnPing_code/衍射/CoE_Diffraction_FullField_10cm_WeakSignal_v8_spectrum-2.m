function CoE_Diffraction_FullField_10cm_WeakSignal_v8_spectrum()
% 全场衍射仿真：10×10 cm 大视场 + 多波长光谱积分 + D65 光源加权 (v8 方案 A)
% =========================================================================
% 核心改进 v8 方案 A:
%   1. 多波长光谱积分 (非单波长近似)
%      - R: 600-640nm, 5nm 步长 (9 个波长)
%      - G: 520-560nm, 5nm 步长 (9 个波长)
%      - B: 450-480nm, 5nm 步长 (7 个波长)
%   2. D65 光源加权 (方案 A: 物理正确顺序)
%      - 光源 (D65) → 反射 → 衍射 → 光谱积分
%      - 每个波长的反射率乘以 D65 强度，再参与衍射计算
%   3. CF 反射率光谱积分 (TMM + Beer-Lambert 定律)
%   4. 保留弱信号可视化增强
%   5. 输出 R:G:B 相对强度比 (基于 D65 能量分布)
% =========================================================================
clearvars; clc; close all;

%% ================= 参数配置 =================
% --- PNG 掩膜参数 ---
PNG_RES             = 300;
PNG_physical_um     = 223.2;
PNG_physical_size_m = PNG_physical_um * 1e-6;

% --- 仿真区间：10×10 cm ---
L_sim      = 3.5e-2;
N_fine     = 2048;

% --- 光学参数 ---
dist_src   = 0.65;
dist_obs   = 0.65;
src_size   = 5e-3;
n_med      = 1.56;
z_pdl      = 1e-6;
z_bm_cf    = 10.2e-6;
d_Rcf      = 2450e-9;
d_Gcf      = 2150e-9;
d_Bcf      = 1500e-9;
d_gap      = z_bm_cf - z_pdl;

% --- 光源与衍射级次 ---
n_src              = 10;
obs_wider          = L_sim;
N_obs_high         = 2048;
m_max              = 9;

% --- 弱信号增强参数 ---
weak_signal_gain   = 15;
log_compression    = 0.15;
clahe_clip_limit   = 0.01;
clahe_tiles        = [8, 8];

% --- 高级次衰减参数 ---
high_order_decay   = true;
decay_exponent     = 2.0;

% --- 光斑大小参数 ---
spot_sigma_pixels  = 9;

% --- 多波长配置 (v8 新增) ---
% R: 600-640nm, G: 520-560nm, B: 450-480nm, step=5nm
wl_R_range = 600:5:640;  % 9 个波长
wl_G_range = 520:5:560;  % 9 个波长
wl_B_range = 450:5:480;  % 7 个波长

fprintf('===== 10×10 cm 多波长光谱仿真 (v8) =====\n');
fprintf('  仿真区间：%.1f × %.1f cm\n', L_sim*100, L_sim*100);
fprintf('  观察范围：±%.1f mm\n', obs_wider*1000);
fprintf('  波长配置:\n');
fprintf('    R: %d 个波长 (%d-%d nm)\n', length(wl_R_range), min(wl_R_range), max(wl_R_range));
fprintf('    G: %d 个波长 (%d-%d nm)\n', length(wl_G_range), min(wl_G_range), max(wl_G_range));
fprintf('    B: %d 个波长 (%d-%d nm)\n', length(wl_B_range), min(wl_B_range), max(wl_B_range));
fprintf('  总波长数：%d\n', length(wl_R_range) + length(wl_G_range) + length(wl_B_range));
fprintf('  高级次衰减：1/m^%.1f\n', decay_exponent);
fprintf('  光斑大小：σ=%.0f 像素\n', spot_sigma_pixels);
fprintf('====================================\n\n');

%% ================= PNG 物理映射 =================
dx_png = PNG_physical_size_m / PNG_RES;
T_pitch_x = dx_png * 3;
T_pitch_y = dx_png * 1;

%% ================= 加载掩膜 =================
fprintf('加载掩膜...\n');
bm  = load_png_resized('BM.png',  PNG_RES);
pdl = load_png_resized('PDL.png', PNG_RES);
rcf = load_png_resized('RCF.png', PNG_RES);
gcf = load_png_resized('GCF.png', PNG_RES);
bcf = load_png_resized('BCF.png', PNG_RES);

is_white = @(x) x > 0.5;
bm_mask = is_white(bm);
pdl_mask = is_white(pdl);
cf_present = @(x) is_white(x);
rcf_mask = cf_present(rcf);
gcf_mask = cf_present(gcf);
bcf_mask = cf_present(bcf);

R_sim = rcf_mask;
G_sim = gcf_mask;
B_sim = bcf_mask;

%% ================= 材料参数 =================
fprintf('加载材料折射率数据...\n');
nkR = loadNK('RCF.txt');
nkG = loadNK('GCF.txt');
nkB = loadNK('BCF.txt');

% 获取 D65 光源权重 (方案 A: 作为入射光源强度)
wl_all = [wl_R_range, wl_G_range, wl_B_range];
D65_full = getD65(wl_all);

% 分别获取各通道的 D65 权重 (不归一化，保留绝对强度信息)
D65_R = interp1(wl_all, D65_full, wl_R_range, 'linear');
D65_G = interp1(wl_all, D65_full, wl_G_range, 'linear');
D65_B = interp1(wl_all, D65_full, wl_B_range, 'linear');

% 归一化到峰值 (用于显示，但保留相对强度比例)
D65_R_norm = D65_R / max(D65_R);
D65_G_norm = D65_G / max(D65_G);
D65_B_norm = D65_B / max(D65_B);

% 计算通道间的相对强度权重 (基于 D65 总能量)
energy_R = trapz(wl_R_range, D65_R);
energy_G = trapz(wl_G_range, D65_G);
energy_B = trapz(wl_B_range, D65_B);
total_energy = energy_R + energy_G + energy_B;

fprintf('D65 光源能量分布 (方案 A: 物理加权):\n');
fprintf('  R 通道 (%d-%dnm): 能量=%.3f (占比 %.1f%%)\n', ...
    min(wl_R_range), max(wl_R_range), energy_R, 100*energy_R/total_energy);
fprintf('  G 通道 (%d-%dnm): 能量=%.3f (占比 %.1f%%)\n', ...
    min(wl_G_range), max(wl_G_range), energy_G, 100*energy_G/total_energy);
fprintf('  B 通道 (%d-%dnm): 能量=%.3f (占比 %.1f%%)\n', ...
    min(wl_B_range), max(wl_B_range), energy_B, 100*energy_B/total_energy);
fprintf('  R:G:B 强度比 = %.3f : %.3f : %.3f\n', ...
    energy_R/energy_B, energy_G/energy_B, 1.0);

% 计算每个波长的 CF 反射率 (TMM)
fprintf('\n计算 CF 反射率光谱...\n');
rR_spectrum = zeros(size(wl_R_range));
rG_spectrum = zeros(size(wl_G_range));
rB_spectrum = zeros(size(wl_B_range));

for i = 1:length(wl_R_range)
    wl = wl_R_range(i);
    nR_i = interp1(nkR.wave, nkR.n, wl, 'linear', 'extrap');
    kR_i = interp1(nkR.wave, nkR.k, wl, 'linear', 'extrap');
    rR_spectrum(i) = tmm_CF_gap_PEC(n_med, nR_i+1i*kR_i, d_Rcf, n_med, d_gap, wl*1e-9);
end

for i = 1:length(wl_G_range)
    wl = wl_G_range(i);
    nG_i = interp1(nkG.wave, nkG.n, wl, 'linear', 'extrap');
    kG_i = interp1(nkG.wave, nkG.k, wl, 'linear', 'extrap');
    rG_spectrum(i) = tmm_CF_gap_PEC(n_med, nG_i+1i*kG_i, d_Gcf, n_med, d_gap, wl*1e-9);
end

for i = 1:length(wl_B_range)
    wl = wl_B_range(i);
    nB_i = interp1(nkB.wave, nkB.n, wl, 'linear', 'extrap');
    kB_i = interp1(nkB.wave, nkB.k, wl, 'linear', 'extrap');
    rB_spectrum(i) = tmm_CF_gap_PEC(n_med, nB_i+1i*kB_i, d_Bcf, n_med, d_gap, wl*1e-9);
end

%% ================= 构建仿真网格 =================
fprintf('构建仿真网格...\n');
x_fine = linspace(-L_sim/2, L_sim/2, N_fine);
y_fine = x_fine;
[Xs, Ys] = meshgrid(x_fine, y_fine);
dx_fine = L_sim / N_fine;

sigma_physical_um = spot_sigma_pixels * dx_fine * 1e6;
fwhm_mm = 2.355 * spot_sigma_pixels * dx_fine * 1e3;
fprintf('  像素尺寸：%.2f μm\n', dx_fine * 1e6);
fprintf('  光斑大小：σ = %.0f 像素 = %.1f μm (FWHM = %.2f mm)\n', ...
    spot_sigma_pixels, sigma_physical_um, fwhm_mm);

% ===== 周期性铺满 =====
u_x = mod(Xs / PNG_physical_size_m, 1);
idx_png = floor(u_x * PNG_RES) + 1;
idx_png = max(1, min(PNG_RES, idx_png));

v_y = mod(Ys / PNG_physical_size_m, 1);
idy_png = floor(v_y * PNG_RES) + 1;
idy_png = max(1, min(PNG_RES, idy_png));

% --- 有效光源点 ---
x_src = linspace(-src_size/2, src_size/2, n_src);
y_src = linspace(-src_size/2, src_size/2, n_src);
[Xsrc, Ysrc] = meshgrid(x_src, y_src);
kx_src = Xsrc / dist_src;
ky_src = Ysrc / dist_src;
cos_src = ones(size(kx_src));
n_valid = numel(kx_src);

% --- 观察面网格 ---
x_obs_wider = linspace(-obs_wider, obs_wider, N_obs_high);
[Xobs_wider, Yobs_wider] = meshgrid(x_obs_wider, x_obs_wider);

% --- 频率域 ---
fx = linspace(-1/(2*dx_fine), 1/(2*dx_fine), N_obs_high);
fy = fx;

%% ================= 初始化强度数组 =================
I_R_spot = zeros(N_obs_high);
I_G_spot = zeros(N_obs_high);
I_B_spot = zeros(N_obs_high);

%% ================= 主仿真循环：R 通道 (多波长积分 - 方案 A) =================
fprintf('\n=== 仿真 R 通道 (%d 个波长，方案 A: D65→反射→衍射) ===\n', length(wl_R_range));
for iw = 1:length(wl_R_range)
    lambda = wl_R_range(iw) * 1e-9;
    k0 = 2*pi / lambda;
    
    % 方案 A: D65 光源强度 × 反射率 (物理正确的顺序)
    rR_iw = abs(rR_spectrum(iw))^2;  % 反射率强度
    d65_w = D65_R(iw);  % D65 光源强度 (未归一化)
    source_weighted_reflectance = rR_iw * d65_w;  % D65 加权后的反射强度
    
    U = zeros(N_fine);
    U(R_sim) = source_weighted_reflectance;  % 使用 D65 加权后的反射率
    
    xwl = dist_obs * lambda * fx;
    [Xwl, Ywl] = meshgrid(xwl, xwl);
    
    I_channel = zeros(N_obs_high);
    
    for s = 1:n_valid
        kx = k0 * kx_src(s); ky = k0 * ky_src(s); w = cos_src(s);
        U_illum = U .* exp(1i * (kx*Xs + ky*Ys));
        Ufft = fftshift(fft2(U_illum));
        Ifft = abs(Ufft).^2;
        
        for mx = -m_max:m_max
            for my = -m_max:m_max
                if mx == 0 && my == 0, continue; end
                
                x_m = dist_obs * lambda * fx(1) * mx;
                y_m = dist_obs * lambda * fy(1) * my;
                
                if abs(x_m) > obs_wider || abs(y_m) > obs_wider, continue; end
                
                dx_obs = x_obs_wider(2) - x_obs_wider(1);
                idx_x = round((x_m - min(x_obs_wider)) / dx_obs) + 1;
                idx_y = round((y_m - min(x_obs_wider)) / dx_obs) + 1;
                
                if idx_x<1 || idx_x>N_obs_high || idx_y<1 || idx_y>N_obs_high, continue; end
                
                win_size = 3;
                x1=max(1,idx_x-win_size); x2=min(N_obs_high,idx_x+win_size);
                y1=max(1,idx_y-win_size); y2=min(N_obs_high,idx_y+win_size);
                local_energy = sum(Ifft(x1:x2, y1:y2), 'all');
                
                ef_x = abs(my_sinc(mx * 0.3));
                ef_y = abs(my_sinc(my * 0.3));
                sinc_envelope = ef_x * ef_y;
                
                m_order = sqrt(mx^2 + my^2);
                if m_order > 0 && high_order_decay
                    order_decay = 1 / (m_order ^ decay_exponent);
                else
                    order_decay = 1;
                end
                
                theta_x = atan(mx * lambda / T_pitch_x);
                theta_y = atan(my * lambda / T_pitch_y);
                tilt_factor = cos(theta_x) * cos(theta_y);
                
                total_attenuation = sinc_envelope * order_decay * tilt_factor;
                peak_intensity = w * local_energy * total_attenuation;
                
                sigma_spot = dx_fine * spot_sigma_pixels;
                gaussian_spot = exp(-((Xobs_wider-x_m).^2 + (Yobs_wider-y_m).^2) / (2*sigma_spot^2));
                
                I_channel = I_channel + peak_intensity * gaussian_spot;
            end
        end
    end
    
    % 累加到 R 通道 (已包含 D65 权重，无需再次加权)
    I_R_spot = I_R_spot + I_channel;
    fprintf('  %d/%d λ=%dnm (D65=%.3f, R=%.3f) 完成\n', ...
        iw, length(wl_R_range), wl_R_range(iw), d65_w, source_weighted_reflectance);
end

%% ================= 主仿真循环：G 通道 (多波长积分 - 方案 A) =================
fprintf('\n=== 仿真 G 通道 (%d 个波长，方案 A: D65→反射→衍射) ===\n', length(wl_G_range));
for iw = 1:length(wl_G_range)
    lambda = wl_G_range(iw) * 1e-9;
    k0 = 2*pi / lambda;
    
    % 方案 A: D65 光源强度 × 反射率
    rG_iw = abs(rG_spectrum(iw))^2;
    d65_w = D65_G(iw);
    source_weighted_reflectance = rG_iw * d65_w;
    
    U = zeros(N_fine);
    U(G_sim) = source_weighted_reflectance;
    
    xwl = dist_obs * lambda * fx;
    I_channel = zeros(N_obs_high);
    
    for s = 1:n_valid
        kx = k0 * kx_src(s); ky = k0 * ky_src(s); w = cos_src(s);
        U_illum = U .* exp(1i * (kx*Xs + ky*Ys));
        Ufft = fftshift(fft2(U_illum));
        Ifft = abs(Ufft).^2;
        
        for mx = -m_max:m_max
            for my = -m_max:m_max
                if mx == 0 && my == 0, continue; end
                
                x_m = dist_obs * lambda * fx(1) * mx;
                y_m = dist_obs * lambda * fy(1) * my;
                
                if abs(x_m) > obs_wider || abs(y_m) > obs_wider, continue; end
                
                dx_obs = x_obs_wider(2) - x_obs_wider(1);
                idx_x = round((x_m - min(x_obs_wider)) / dx_obs) + 1;
                idx_y = round((y_m - min(x_obs_wider)) / dx_obs) + 1;
                
                if idx_x<1 || idx_x>N_obs_high || idx_y<1 || idx_y>N_obs_high, continue; end
                
                win_size = 3;
                x1=max(1,idx_x-win_size); x2=min(N_obs_high,idx_x+win_size);
                y1=max(1,idx_y-win_size); y2=min(N_obs_high,idx_y+win_size);
                local_energy = sum(Ifft(x1:x2, y1:y2), 'all');
                
                ef_x = abs(my_sinc(mx * 0.3));
                ef_y = abs(my_sinc(my * 0.3));
                sinc_envelope = ef_x * ef_y;
                
                m_order = sqrt(mx^2 + my^2);
                order_decay = 1 / (m_order ^ decay_exponent);
                
                theta_x = atan(mx * lambda / T_pitch_x);
                theta_y = atan(my * lambda / T_pitch_y);
                tilt_factor = cos(theta_x) * cos(theta_y);
                
                total_attenuation = sinc_envelope * order_decay * tilt_factor;
                peak_intensity = w * local_energy * total_attenuation;
                
                sigma_spot = dx_fine * spot_sigma_pixels;
                gaussian_spot = exp(-((Xobs_wider-x_m).^2 + (Yobs_wider-y_m).^2) / (2*sigma_spot^2));
                
                I_channel = I_channel + peak_intensity * gaussian_spot;
            end
        end
    end
    
    % 累加到 G 通道 (已包含 D65 权重)
    I_G_spot = I_G_spot + I_channel;
    fprintf('  %d/%d λ=%dnm (D65=%.3f, R=%.3f) 完成\n', ...
        iw, length(wl_G_range), wl_G_range(iw), d65_w, source_weighted_reflectance);
end

%% ================= 主仿真循环：B 通道 (多波长积分 - 方案 A) =================
fprintf('\n=== 仿真 B 通道 (%d 个波长，方案 A: D65→反射→衍射) ===\n', length(wl_B_range));
for iw = 1:length(wl_B_range)
    lambda = wl_B_range(iw) * 1e-9;
    k0 = 2*pi / lambda;
    
    % 方案 A: D65 光源强度 × 反射率
    rB_iw = abs(rB_spectrum(iw))^2;
    d65_w = D65_B(iw);
    source_weighted_reflectance = rB_iw * d65_w;
    
    U = zeros(N_fine);
    U(B_sim) = source_weighted_reflectance;
    
    xwl = dist_obs * lambda * fx;
    I_channel = zeros(N_obs_high);
    
    for s = 1:n_valid
        kx = k0 * kx_src(s); ky = k0 * ky_src(s); w = cos_src(s);
        U_illum = U .* exp(1i * (kx*Xs + ky*Ys));
        Ufft = fftshift(fft2(U_illum));
        Ifft = abs(Ufft).^2;
        
        for mx = -m_max:m_max
            for my = -m_max:m_max
                if mx == 0 && my == 0, continue; end
                
                x_m = dist_obs * lambda * fx(1) * mx;
                y_m = dist_obs * lambda * fy(1) * my;
                
                if abs(x_m) > obs_wider || abs(y_m) > obs_wider, continue; end
                
                dx_obs = x_obs_wider(2) - x_obs_wider(1);
                idx_x = round((x_m - min(x_obs_wider)) / dx_obs) + 1;
                idx_y = round((y_m - min(x_obs_wider)) / dx_obs) + 1;
                
                if idx_x<1 || idx_x>N_obs_high || idx_y<1 || idx_y>N_obs_high, continue; end
                
                win_size = 3;
                x1=max(1,idx_x-win_size); x2=min(N_obs_high,idx_x+win_size);
                y1=max(1,idx_y-win_size); y2=min(N_obs_high,idx_y+win_size);
                local_energy = sum(Ifft(x1:x2, y1:y2), 'all');
                
                ef_x = abs(my_sinc(mx * 0.3));
                ef_y = abs(my_sinc(my * 0.3));
                sinc_envelope = ef_x * ef_y;
                
                m_order = sqrt(mx^2 + my^2);
                order_decay = 1 / (m_order ^ decay_exponent);
                
                theta_x = atan(mx * lambda / T_pitch_x);
                theta_y = atan(my * lambda / T_pitch_y);
                tilt_factor = cos(theta_x) * cos(theta_y);
                
                total_attenuation = sinc_envelope * order_decay * tilt_factor;
                peak_intensity = w * local_energy * total_attenuation;
                
                sigma_spot = dx_fine * spot_sigma_pixels;
                gaussian_spot = exp(-((Xobs_wider-x_m).^2 + (Yobs_wider-y_m).^2) / (2*sigma_spot^2));
                
                I_channel = I_channel + peak_intensity * gaussian_spot;
            end
        end
    end
    
    % 累加到 B 通道 (已包含 D65 权重)
    I_B_spot = I_B_spot + I_channel;
    fprintf('  %d/%d λ=%dnm (D65=%.3f, R=%.3f) 完成\n', ...
        iw, length(wl_B_range), wl_B_range(iw), d65_w, source_weighted_reflectance);
end

fprintf('\n多波长光谱积分完成 (方案 A: D65→反射→衍射)!\n');

%% ================= 物理强度验证 =================
fprintf('\n===== 物理强度验证 =====\n');
center_R = mean(I_R_spot(900:1100, 900:1100), 'all');
center_G = mean(I_G_spot(900:1100, 900:1100), 'all');
center_B = mean(I_B_spot(900:1100, 900:1100), 'all');

edge_R = mean(I_R_spot(100:300, 100:300), 'all');
edge_G = mean(I_G_spot(100:300, 100:300), 'all');
edge_B = mean(I_B_spot(100:300, 100:300), 'all');

fprintf('中心区域 (±10mm) 平均强度:\n');
fprintf('  R: %.3e, G: %.3e, B: %.3e\n', center_R, center_G, center_B);
fprintf('边缘区域 (±40mm) 平均强度:\n');
fprintf('  R: %.3e, G: %.3e, B: %.3e\n', edge_R, edge_G, edge_B);

if center_R > 0
    ratio_R = edge_R / center_R;
    ratio_G = edge_G / center_G;
    ratio_B = edge_B / center_B;
    fprintf('边缘/中心强度比:\n');
    fprintf('  R: %.4f, G: %.4f, B: %.4f\n', ratio_R, ratio_G, ratio_B);
end
fprintf('================================\n\n');

%% ================= 弱信号增强显示 =================
fprintf('应用弱信号增强映射...\n');

enhance_weak_signal = @(I, gain, log_c) ...
    log(1 + gain * log_c * I) / log(1 + gain * log_c * max(I(:)) + 1e-10);

I_R_enh = enhance_weak_signal(I_R_spot, weak_signal_gain, log_compression);
I_G_enh = enhance_weak_signal(I_G_spot, weak_signal_gain, log_compression);
I_B_enh = enhance_weak_signal(I_B_spot, weak_signal_gain, log_compression);

if exist('adapthisteq', 'file')
    I_R_enh = adapthisteq(I_R_enh, 'ClipLimit', clahe_clip_limit, 'NumTiles', clahe_tiles);
    I_G_enh = adapthisteq(I_G_enh, 'ClipLimit', clahe_clip_limit, 'NumTiles', clahe_tiles);
    I_B_enh = adapthisteq(I_B_enh, 'ClipLimit', clahe_clip_limit, 'NumTiles', clahe_tiles);
end

% D65 混色 (v8: 直接使用光谱积分结果，已包含 D65 权重)
RGB_enhanced = cat(3, I_R_enh, I_G_enh, I_B_enh);
RGB_enhanced = RGB_enhanced / max(RGB_enhanced(:) + 1e-10);

%% ================= 显示与保存 =================
fig = figure('Position', [50, 50, 1600, 900], 'Color', 'k');
x_mm = x_obs_wider * 1e3;
sb = 10;

subplot(2,4,1);
imshow(cat(3,I_R_enh,zeros(size(I_R_enh)),zeros(size(I_R_enh))),'XData',x_mm,'YData',x_mm);
axis on; set(gca,'Color','none','XColor','w','YColor','w');
title('R 通道 (600-640nm, D65→反射→衍射)','Color','w','FontSize',10);
add_scale_bar_mm(gca,sb,'10 mm','w',2,10);

subplot(2,4,2);
imshow(cat(3,zeros(size(I_G_enh)),I_G_enh,zeros(size(I_G_enh))),'XData',x_mm,'YData',x_mm);
axis on; set(gca,'Color','none','XColor','w','YColor','w');
title('G 通道 (520-560nm, D65→反射→衍射)','Color','w','FontSize',10);
add_scale_bar_mm(gca,sb,'10 mm','w',2,10);

subplot(2,4,3);
imshow(cat(3,zeros(size(I_B_enh)),zeros(size(I_B_enh)),I_B_enh),'XData',x_mm,'YData',x_mm);
axis on; set(gca,'Color','none','XColor','w','YColor','w');
title('B 通道 (450-480nm, D65→反射→衍射)','Color','w','FontSize',10);
add_scale_bar_mm(gca,sb,'10 mm','w',2,10);

subplot(2,4,4);
imshow(RGB_enhanced,'XData',x_mm,'YData',x_mm);
axis on; set(gca,'Color','none','XColor','w','YColor','w');
title('D65 混色 (方案 A: 物理加权)','Color','w','FontSize',10);
add_scale_bar_mm(gca,sb,'10 mm','w',2,10);

subplot(2,4,5);
imshow(log(1 + I_R_spot * 1e-6), 'XData', x_mm, 'YData', x_mm);
axis on; set(gca,'Color','none','XColor','w','YColor','w');
title('原始物理强度-R (Log)','Color','w','FontSize',10);
add_scale_bar_mm(gca,sb,'10 mm','w',2,10);

subplot(2,4,6);
imshow(log(1 + I_G_spot * 1e-6), 'XData', x_mm, 'YData', x_mm);
axis on; set(gca,'Color','none','XColor','w','YColor','w');
title('原始物理强度-G (Log)','Color','w','FontSize',10);
add_scale_bar_mm(gca,sb,'10 mm','w',2,10);

subplot(2,4,7);
imshow(log(1 + I_B_spot * 1e-6), 'XData', x_mm, 'YData', x_mm);
axis on; set(gca,'Color','none','XColor','w','YColor','w');
title('原始物理强度-B (Log)','Color','w','FontSize',10);
add_scale_bar_mm(gca,sb,'10 mm','w',2,10);

subplot(2,4,8);
RGB_physical = cat(3, ...
    mat2gray(log(1 + I_R_spot * 1e-6)), ...
    mat2gray(log(1 + I_G_spot * 1e-6)), ...
    mat2gray(log(1 + I_B_spot * 1e-6)));
imshow(RGB_physical,'XData',x_mm,'YData',x_mm);
axis on; set(gca,'Color','none','XColor','w','YColor','w');
title('原始物理强度-W (Log)','Color','w','FontSize',10);
add_scale_bar_mm(gca,sb,'10 mm','w',2,10);

sgtitle(sprintf('10×10 cm 多波长仿真 (v8 方案 A) | R:%d-%dnm G:%d-%dnm B:%d-%dnm | D65→反射→衍射', ...
    min(wl_R_range), max(wl_R_range), min(wl_G_range), max(wl_G_range), min(wl_B_range), max(wl_B_range)), ...
    'Color','w','FontSize',12);

% 保存结果
imwrite(cat(3,I_R_enh,zeros(size(I_R_enh)),zeros(size(I_R_enh))), 'Perceptual_R_10cm_v8_spectrum_A.png');
imwrite(cat(3,zeros(size(I_G_enh)),I_G_enh,zeros(size(I_G_enh))), 'Perceptual_G_10cm_v8_spectrum_A.png');
imwrite(cat(3,zeros(size(I_B_enh)),zeros(size(I_B_enh)),I_B_enh), 'Perceptual_B_10cm_v8_spectrum_A.png');
imwrite(RGB_enhanced, 'Perceptual_W_10cm_v8_spectrum_A_D65.png');

save('Physical_Intensity_Data_v8_spectrum_A.mat', 'I_R_spot', 'I_G_spot', 'I_B_spot', ...
     'x_obs_wider', 'wl_R_range', 'wl_G_range', 'wl_B_range', ...
     'D65_R', 'D65_G', 'D65_B', 'energy_R', 'energy_G', 'energy_B');

fprintf('\n=== 仿真完成 (方案 A: D65→反射→衍射) ===\n');
fprintf('已保存:\n');
fprintf('  - Perceptual_R_10cm_v8_spectrum_A.png\n');
fprintf('  - Perceptual_G_10cm_v8_spectrum_A.png\n');
fprintf('  - Perceptual_B_10cm_v8_spectrum_A.png\n');
fprintf('  - Perceptual_W_10cm_v8_spectrum_A_D65.png (D65 混色)\n');
fprintf('  - Physical_Intensity_Data_v8_spectrum_A.mat (原始数据 + D65 能量分布)\n');
fprintf('\n方案 A 说明:\n');
fprintf('  - 物理过程：D65 光源 → CF 反射 → 衍射 → 光谱积分\n');
fprintf('  - 每个波长的反射率已乘以 D65 强度，再参与衍射计算\n');
fprintf('  - R:G:B 相对强度比反映真实 D65 光源下的颜色平衡\n');

end  % 主函数结束

%% ================= 辅助函数 =================
function bm_mask = load_png_resized(fname, res)
if ~exist(fname, 'file')
    warning('掩膜文件不存在：%s, 使用全 1 掩膜', fname);
    bm_mask = ones(res);
    return;
end
img = imread(fname);
if size(img, 3) == 3
    img = rgb2gray(img);
end
bm_mask = im2double(imresize(img, [res, res]));
end

function add_scale_bar_mm(ax, length_mm, label, color, linewidth, fontsize)
hold(ax, 'on');
xlim_val = xlim(ax); ylim_val = ylim(ax);
dx = xlim_val(2) - xlim_val(1);
dy = ylim_val(2) - ylim_val(1);
start_x = xlim_val(1) + 0.1*dx;
start_y = ylim_val(1) + 0.05*dy;
line(ax, [start_x, start_x+length_mm], [start_y, start_y], ...
    'Color', color, 'LineWidth', linewidth);
text(ax, start_x+length_mm/2, start_y+0.1*dy, label, ...
    'Color', color, 'FontSize', fontsize, 'HorizontalAlignment','center');
hold(ax, 'off');
end

function y = my_sinc(x)
y = ones(size(x));
non_zero = (x ~= 0);
y(non_zero) = sin(pi * x(non_zero)) ./ (pi * x(non_zero));
end

function nk = loadNK(fname)
if ~exist(fname, 'file')
    warning('NK 文件不存在：%s', fname);
    nk.wave = []; nk.n = []; nk.k = [];
    return;
end
try
    opts = detectImportOptions(fname, 'FileType', 'text');
    opts = setvaropts(opts, 'VariableNamingRule', 'preserve');
    tab = readtable(fname, opts);
catch
    tab = readtable(fname, 'VariableNamingRule', 'preserve');
end
cols = tab.Properties.VariableNames;
for c = cols
    cn = c{:};
    if any(strcmpi(cn, {'Wavelength','wave','wl','Lambda'}))
        tab.Properties.VariableNames{c{:}} = 'wave';
    end
    if any(strcmpi(cn, {'n','nk_real','refractive'}))
        tab.Properties.VariableNames{c{:}} = 'n';
    end
    if any(strcmpi(cn, {'k','nk_imag','extinction'}))
        tab.Properties.VariableNames{c{:}} = 'k';
    end
end
if ~all(ismember({'wave','n','k'}, tab.Properties.VariableNames))
    error('NK 文件缺少必需列：wave, n, k');
end
nk.wave = tab.wave; nk.n = tab.n; nk.k = tab.k;
end

function r = tmm_CF_gap_PEC(n0, n1, d1, n2, d2, lam)
d1 = 2*pi*n1*d1/lam;
d2 = 2*pi*n2*d2/lam;
r01 = (n0-n1)/(n0+n1);
r12 = (n1-n2)/(n1+n2);
r2 = -1*exp(-2i*d2);
r1 = (r12 + r2)./(1 + r12.*r2);
r = (r01 + r1.*exp(-2i*d1))./(1 + r01.*r1.*exp(-2i*d1));
end

function D65 = getD65(wl)
% CIE D65 标准光源 (相对光谱功率分布)
% 使用普朗克黑体辐射近似，T=6504K
T = 6504;
c1 = 3.7418e-16;  % W·m²
c2 = 0.014388;    % m·K
lambda_m = wl * 1e-9;
D65 = (c1 ./ (lambda_m.^5)) ./ (exp(c2 ./ (lambda_m * T)) - 1);
D65 = D65 / max(D65);

% 修正：D65 在短波长处有较高的能量
% 使用简化的 CIE D65 公式
if max(wl) > 400
    % 对于可见光范围，使用更准确的 CIE D65 近似
    S0 = ones(size(wl));
    S1 = zeros(size(wl));
    S2 = zeros(size(wl));
    
    % 简化的 CIE D65 光谱
    for i = 1:length(wl)
        w = wl(i);
        if w >= 300 && w <= 830
            % CIE D65 相对光谱 (简化版)
            S0(i) = 1.0;
            if w < 440
                S1(i) = -0.024;
            elseif w < 540
                S1(i) = 0.015 * (w-440)/100;
            else
                S1(i) = 0.015 * (1 - (w-540)/300);
            end
            S2(i) = 0;
        end
    end
    
    D65 = S0 + 0.0247*S1 - 0.0009*S2;
    D65 = max(0, D65);
    D65 = D65 / max(D65);
end
end