
% function  [] = detectIEDs_overDirs(dirName)
dirName = '/mnt/mfs/patients/cereplex/CUCX2/';

% finding all of the files
if strcmp(dirName(end),'/')
    files = subdir([dirName '*.ns5']);
else
    files = subdir([dirName '/*.ns5']);
end

% looping over files
for fl = 1:length(files)
    detectIEDs_fromFile(files(fl).name)
end
