clear, clc, close all

lambda = 550e-9;
z = 0.7;
Pitch = 2e-3;
r = 0.5e-3;
nn = 700;
Nx_period = 2;
Ny_period = 2;
Fx = 5e-3;
Fy = 5e-3;

x_cell = linspace(-Pitch/2, Pitch/2, nn);
y_cell = linspace(-Pitch/2, Pitch/2, nn);
[Xc, Yc] = meshgrid(x_cell, y_cell);

xmin = -Nx_period * Pitch / 2;
xmax = -xmin;
ymin = -Ny_period * Pitch / 2;
ymax = -ymin;
Xmin = -Fx;
Xmax = Fx;
Ymin = -Fy;
Ymax = Fy;
Tnn = Nx_period * nn;

x = linspace(xmin, xmax, Tnn);
y = linspace(ymin, ymax, Tnn);
[x, y] = meshgrid(x, y);
X = linspace(Xmin, Xmax, Tnn);
Y = linspace(Ymin, Ymax, Tnn);
[X, Y] = meshgrid(X, Y);

cellunit = 2;
SinglePattern1 = CustomPattern(0, Xc, Yc, r, r, 0, 0);
SinglePattern2 = CustomPattern(0, Xc, Yc, r, r, 0, 0);
SinglePattern3 = CustomPattern(0, Xc, Yc, r, r, 0, 0);
SinglePattern4 = CustomPattern(0, Xc, Yc, r, r, 0, 0);
cellPattern = [SinglePattern1 SinglePattern2; SinglePattern3 SinglePattern4];
Mask = repmat(cellPattern, Ny_period/cellunit, Nx_period/cellunit);

k = 2 * pi / lambda;
[~, U2] = fraunhofer_fft_SPH(Mask, xmin, xmax, ymin, ymax, Tnn, Tnn, ...
    lambda, z, Xmin, Xmax, Ymin, Ymax, Tnn, Tnn, k, 1);

I = abs(U2).^2;
I_norm = I ./ max(I(:));
I_sorted = sort(I_norm(:));
clipLevel = I_sorted(max(1, round(0.9985 * numel(I_sorted))));
I_clip = min(I_norm, clipLevel) ./ clipLevel;
alpha = 60;
gamma = 0.72;
I_vis = log(1 + alpha * I_clip) / log(1 + alpha);
I_vis = I_vis .^ gamma;

outDir = fullfile('..', 'ImageForShow');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

fig = figure('Color', 'w', 'Position', [100 100 760 1100]);
subplot(2, 1, 1)
h1 = imagesc(x(1,:), y(:,1), Mask);
set(h1, 'Interpolation', 'nearest');
axis image
set(gca, 'YDir', 'normal');
colormap gray
xlabel('x (m)')
ylabel('y (m)')
title('Four circular apertures')

subplot(2, 1, 2)
imagesc(X(1,:)*1e3, Y(:,1)*1e3, I_vis);
axis image
set(gca, 'YDir', 'normal');
colormap gray
colorbar
clim([0 1])
xlabel('x / mm')
ylabel('y / mm')
title('Four-aperture diffraction')

exportgraphics(fig, fullfile(outDir, 'four_aperture_preview.png'), 'Resolution', 200);
save(fullfile(outDir, 'four_aperture_preview_data.mat'), ...
    'Mask', 'I', 'I_norm', 'I_vis', 'x', 'y', 'X', 'Y', ...
    'lambda', 'z', 'Pitch', 'r', 'Fx', 'Fy', 'nn', 'clipLevel', 'alpha', 'gamma');
