function export_matlab_colormaps()
%EXPORT_MATLAB_COLORMAPS Export MATLAB's standard maps for offline viewing.
    names = ["parula","turbo","hsv","hot","cool","spring","summer","autumn", ...
             "winter","gray","bone","copper","pink","sky","abyss","jet", ...
             "lines","colorcube","prism","flag","white"];
    if exist('colormaplist','file') || exist('colormaplist','builtin')
        names = unique([names, reshape(string(colormaplist),1,[])],'stable');
    end
    for name = names
        map = feval(char(name),256);
        assert(isequal(size(map),[256,3]),'Unexpected map size for %s',name);
        palettes.(name) = uint8(round(255 * map));
    end
    payload = struct('release',version('-release'),'palettes',palettes);
    path = fullfile(fileparts(mfilename('fullpath')),'matlab_colormaps.json');
    fid = fopen(path,'w');
    assert(fid >= 0,'Could not open %s',path);
    closeFile = onCleanup(@() fclose(fid));
    fprintf(fid,'%s\n',jsonencode(payload));
end
