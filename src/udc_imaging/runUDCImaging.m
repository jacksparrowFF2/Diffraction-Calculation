function [result, images] = runUDCImaging(cfg, scene, encoding)
%runUDCImaging Pupil-plane transmission, single-pass phase and image formation.
% scene: image-plane irradiance at cfg.sensorPitch (not an object-plane image).
% encoding="linear" accepts HDR; "srgb" explicitly decodes an ordinary image.
    arguments
        cfg (1,1) struct = udcImagingDefaults()
        scene {mustBeNumericOrLogical} = []
        encoding (1,1) string {mustBeMember(encoding, ["linear", "srgb"])} = "linear"
    end
    assert(cfg.mode == "udc", 'DisplayDiffraction:Mode', 'Use udc configuration.');
    if isempty(scene)
        moduleDir = fileparts(mfilename('fullpath'));
        scene = imread(fullfile(moduleDir, 'input', 'USAF500.png'));
        encoding = "srgb";
    end
    result = computeUDCOptics(cfg);
    images = formUDCImage(scene, result, encoding);
end
