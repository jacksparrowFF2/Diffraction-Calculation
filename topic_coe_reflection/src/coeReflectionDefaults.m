function cfg = coeReflectionDefaults()
%coeReflectionDefaults Create a complete COE reflection configuration in SI.
    arguments
    end
    cfg = displayDiffractionDefaults("reflection");
    cfg.outputSize = 1201;
end
