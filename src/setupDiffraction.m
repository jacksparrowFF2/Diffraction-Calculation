function projectRoot = setupDiffraction()
%setupDiffraction Add public feature folders to MATLAB's search path.
% Call after addpath('src'). MATLAB automatically resolves each source folder's private dependencies.
% Imported COE reference programs are opt-in; see topic_coe_reference/docs/README.md.
    arguments
    end
    sourceDir = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(sourceDir);
    modules = ["topic_aperture_fraunhofer", "topic_phase_layer", "topic_coe_reflection", ...
        "topic_udc_imaging", "topic_pixel_unfolding"];
    addpath(sourceDir);
    for module = modules
        addpath(fullfile(projectRoot, module, "src"));
    end
end
