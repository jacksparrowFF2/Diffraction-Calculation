function cfg = displayDiffractionDefaults(mode)
%displayDiffractionDefaults Compatibility entry for feature-specific defaults.
% Prefer coeReflectionDefaults or udcImagingDefaults for standalone feature use.
    arguments
        mode (1,1) string {mustBeMember(mode, ["reflection", "udc"])} = "reflection"
    end
    previousPath = path;
    restorePath = onCleanup(@() path(previousPath));
    setupDiffraction();
    if mode == "reflection"
        cfg = coeReflectionDefaults();
    else
        cfg = udcImagingDefaults();
    end
end
