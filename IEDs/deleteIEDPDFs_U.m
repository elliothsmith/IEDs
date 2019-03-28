function [] = deleteIEDPDFs_U()
% DELETEIEDPDFS deletes pdfs generated from IEDwaves_U in order to
%   re-generate those pdfs.

% author: EHS20181120

% where to save figs.
figDir = '/media/user1/data4TB/Figs/IEDs';

% looping over patients
nPts = 5;
for pt = 1:6
    % (s)which patients
    switch pt
        case 1
            ptID = 'nu';
            % sad magic variables resulting from my lack of prescience go here
            Fs = 3e4;
        case 2
            ptID = 'gamma';
            Fs = 3e4;
        case 3
            ptID = 'epsilon';
            Fs = 3e4;
        case 4
            ptID = 'eta';
            Fs = 3e4;
        case 5
            ptID = 'kappa';
            Fs = 3e4;
        case 6
            ptID = 'zeta';
            Fs = 3e4;
    end
    
    
    % getting lists of files
    dirList = dir(fullfile(figDir,[ptID '/*.pdf']));
    
    for fl = 1:length(dirList)
                
        delete(fullfile(dirList(fl).folder,dirList(fl).name))
        
    end
end