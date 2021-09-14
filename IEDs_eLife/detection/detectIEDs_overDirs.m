function  [] = detectIEDs_overDirs(ptID,dirName);

% finding all of the files
if strcmp(dirName(end),'/')
    files = subdir([dirName ptID '/*.ns5']);
else
    files = subdir([dirName '/' ptID '/*.ns5']);
end
    
% looping over files
for fl = startFile:length(files)
    fprintf('\nfile: %s',files(fl).name)
    detectIEDs_fromFile(ptID,files(fl).name)
end


