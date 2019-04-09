% function  [] = detectIEDs_overDirs(dirName)
nPts = 4;
dirName = '/mnt/mfs/patients/cereplex/';
for pts = 4
    % looping over patients and defining the startFile in case some files have already been processed. 
    switch pts
        case 1
            ptID = 'CUCX2';
			startFile = 1;
        case 2
            ptID = 'CUCX3';
			startFile = 188;
			skipFiles = [93 106 186 187];
        case 3
            ptID = 'CUCX4';
			startFile = 1;
		case 4 
            ptID = 'CUCX5';
			startFile = 1;
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
