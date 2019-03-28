% function  [] = detectIEDs_overDirs(dirName)
nPts = 13;
try
    dirName = '/home/user1/Desktop/R60/Data/UEA/UofU/';
catch
    % sudo mount.cifs //ClinicalWorkstation/ClinicalWorkstation_R60 /home/user1/Desktop/R60 -o user=Administrator
end  % I realize this doesn't make sense, but I don't care!!
for pts = 1:13
    % looping over patients
    switch pts
        case 1
            ptID = 'nu';
        case 2
            ptID = 'gamma';
        case 3
            ptID = 'beta';
        case 4
            ptID = 'delta';
        case 5
            ptID = 'epsilon';
        case 6
            ptID = 'eta';
        case 7
            ptID = 'iota';
        case 8
            ptID = 'kappa';
        case 9
            ptID = 'lambda';
        case 10
            ptID = 'mu';
        case 11
            ptID = 'theta';
        case 12
            ptID = 'xi';
        case 13
            ptID = 'zeta';
    end
    
    % finding all of the files
    if strcmp(dirName(end),'/')
        files = subdir([dirName ptID '/*.ns5']);
    else
        files = subdir([dirName '/' ptID '/*.ns5']);
    end
    
    % looping over files
    for fl = 1:length(files)
        fprintf('\nfile: %s',files(fl).name)
        detectIEDs_fromFile(ptID,files(fl).name)
    end
    
end
