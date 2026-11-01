function projectRoot = setupDiffraction()
%setupDiffraction Add public feature folders to MATLAB's search path.
% Call after addpath('src'). Private dependencies remain scoped to each feature.
% Imported COE reference programs are opt-in; see coe_reference/README.md.
    arguments
    end
    sourceDir = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(sourceDir);
    modules = ["aperture_diffraction", "phase_layer", "coe_reflection", ...
        "udc_imaging", "pixel_unfolding"];
    addpath(sourceDir);
    for module = modules
        addpath(fullfile(sourceDir, module));
    end
end
