
function moveRawFiles


dirPath = uigetdir;
rawFolder = createFolder([dirPath '\raw\']);
movefile(fullfile(dirPath, 'file*.mat'), rawFolder);





