function outputFolder = demoCOEUDC(outputFolder)
%demoCOEUDC Run both branches, export previews, physical arrays and metrics.
% All optical coefficients are example values. Plots use a shared display scale.
    arguments
        outputFolder (1,1) string = string(fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
            'ImageForShow', 'COE_UDC'))
    end
    if ~isfolder(outputFolder)
        mkdir(outputFolder);
    end
    reflectionCfg = displayDiffractionDefaults("reflection");
    % Reproducible synthetic height map. Replace with h_map from the phase-layer
    % script or a measured map of exactly the same grid and physical extent.
    n = reflectionCfg.periods*reflectionCfg.samplesPerPeriod;
    [x, y] = meshgrid((0:n(1)-1)*reflectionCfg.pitch/reflectionCfg.samplesPerPeriod, ...
        (0:n(2)-1)*reflectionCfg.pitch/reflectionCfg.samplesPerPeriod);
    reflectionCfg.heightMap = 110e-9*(1 + sin(2*pi*x/17.3e-6).*sin(2*pi*y/23.7e-6));
    reflection = runCOEReflection(reflectionCfg);
    udcCfg = displayDiffractionDefaults("udc");
    udcCfg.heightMap = reflectionCfg.heightMap;
    [udc, images] = runUDCImaging(udcCfg);
    writetable(reflection.metrics, fullfile(outputFolder, 'reflection_metrics.csv'));
    writetable(udc.metrics, fullfile(outputFolder, 'udc_metrics.csv'));

    fig = figure('Name', 'COE reflection and UDC', 'Color', 'w', 'Position', [80, 80, 1400, 850]);
    tiledlayout(fig, 2, 3, 'TileSpacing', 'compact');
    nexttile;
    imagesc(reflection.x*1e6, reflection.y*1e6, reflection.mask);
    axis image;
    title('Shared aperture (example)');
    xlabel('x / um'); ylabel('y / um');
    nexttile;
    referencePeak = max(reflection.baselineKernel(:));
    imagesc(rad2deg(reflection.coordinate), rad2deg(reflection.coordinate), ...
        log10(max(reflection.kernel(:, :, 2)/referencePeak, 1e-7)));
    axis image; clim([-7, 0]); colorbar;
    title('Reflection, green: log10(I / baseline peak)');
    xlabel('Paraxial angle / deg'); ylabel('Paraxial angle / deg');
    nexttile;
    centre = (udcCfg.outputSize + 1)/2;
    plot(udc.spatialFrequency(centre:end), squeeze(udc.mtf(centre, centre:end, :)), 'LineWidth', 1.2);
    ylim([0, 1]); grid on;
    xlabel('Spatial frequency / cycles per mm'); ylabel('MTF');
    title('UDC horizontal MTF (finite PSF, pixel integrated)');
    legend('R proxy', 'G proxy', 'B proxy');
    nexttile;
    image(encodeSrgb(images.clear)); axis image off;
    title('Clear pupil, same exposure');
    nexttile;
    image(encodeSrgb(images.fixedExposure)); axis image off;
    title('UDC, same exposure');
    nexttile;
    image(encodeSrgb(images.brightnessMatched)); axis image off;
    title('UDC, brightness matched (optics only)');
    exportgraphics(fig, fullfile(outputFolder, 'coe_udc_overview.png'), 'Resolution', 160);
    save(fullfile(outputFolder, 'coe_udc_results.mat'), 'reflection', 'udc', 'images', '-v7.3');
    disp(reflection.metrics);
    disp(udc.metrics);
    fprintf('Outputs: %s\n', outputFolder);
end

function encoded = encodeSrgb(linear)
    linear = min(max(linear, 0), 1); % Display ONLY; physical arrays retain HDR.
    encoded = 12.92*linear;
    high = linear > 0.0031308;
    encoded(high) = 1.055*linear(high).^(1/2.4) - 0.055;
end
