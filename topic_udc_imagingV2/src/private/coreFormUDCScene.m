function images = coreFormUDCScene(sceneImage,optics,encoding)
%coreFormUDCScene 使用 V1 的线性传感器成像函数。
    srcDir = fileparts(fileparts(mfilename('fullpath')));
    repositoryDir = fileparts(fileparts(srcDir));
    savedPath = path;
    restorePath = onCleanup(@() path(savedPath));
    addpath(fullfile(repositoryDir,'topic_udc_imaging','src'));
    outputSize = [];
    frame = struct();
    if isfield(optics.config,'sceneProjection')
        projection = optics.config.sceneProjection;
        outputSize = fliplr(projection.imageSizePixels);
        if isfield(projection,'frame'), frame = projection.frame; end
    end
    images = formUDCImage(sceneImage,optics,encoding,outputSize,frame);
end
