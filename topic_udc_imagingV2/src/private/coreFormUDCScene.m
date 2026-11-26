function images = coreFormUDCScene(sceneImage,optics,encoding)
%coreFormUDCScene 使用 V1 的线性传感器成像函数。
    srcDir = fileparts(fileparts(mfilename('fullpath')));
    repositoryDir = fileparts(fileparts(srcDir));
    savedPath = path;
    restorePath = onCleanup(@() path(savedPath));
    addpath(fullfile(repositoryDir,'topic_udc_imaging','src'));
    images = formUDCImage(sceneImage,optics,encoding);
end
