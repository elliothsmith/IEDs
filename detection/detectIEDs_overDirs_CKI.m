% function  [] = detectIEDs_overDirs(dirName)
nPts = 5;
dirName = '/mnt/mfs/patients/CKI/';
for pts = 4:nPts % c2 and c4 don't appear to be mounted...
	clearvars -except nPts dirName pts 
    % looping over patients
    switch pts
        case 1
            ptID = 'c2';
        case 2
            ptID = 'c3';
        case 3
            ptID = 'c4';
        case 4
            ptID = 'c5';
			startFile = 10;
        case 5
            ptID = 'c7';
			startFile = 900;
    end
    
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
    
end
