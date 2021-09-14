% [20190520] figuring out dates in order to plot IEDs over time.
% starting with C7 - need to loop over patients.
% [20190718] update: put in patient loop and seizure data
%	- seizure onsets and offsets aree in samples.
% ^^^  from previous iterations implemented on saturn for CUMC patients. ^^^
% [20190827] now getting this running on the utah data.

clear;

% as long as I'm hard coding this, I may as well put it right up front..
Fs = 3e4;

close all;

% looping over patients.
for pt = 4
    clearvars -except Fs pt circStats HR_p HR_t R_z R_p k_concParam
    switch pt
        case 1
            patientID = 'c3';
        case 2
            patientID = 'c4';
		case 3
			patientID = 'c5';
		case 4
			patientID = 'c7';
		case 5
			patientID = 'CUCX2';
		case 6
			patientID = 'CUCX3';
		case 7
			patientID = 'CUCX4';
    end
	nIEDs = 0;
   
    % for the Utah data we need to loop over the files in the directory.
    dirList = dir(sprintf('/storage/selected_data/elliotWorking/Data/IEDs/allIEDs_%s/*.mat',patientID));
    for fl = 1:length(dirList)
        load(fullfile(dirList(fl).folder,dirList(fl).name))
        
        % loop over files.
        nIEDsFl = length(thisIED);
        for id = 1:nIEDsFl
            updateUser('compiling IED timing data for IED',nIEDs + id,100,nIEDs)
           
            % other features in the data we may want to visualaize.
            sigIEDs(nIEDs + id) = thisIED(id).wavePval<0.05;
            V = pinv(thisIED(id).waveBetas(1:2));
            IEDspeed(nIEDs + id) = sqrt(sum(V.^2));
            IEDdirection(nIEDs + id) = angle(V(1) + V(2)*sqrt(-1));
            clear V

			IEDwf(nIEDs + id,:) = thisIED(id).medianWF;
            
        end % looping over IEDs in each file.
        nIEDs = nIEDs + nIEDsFl;
    end % looping over files.  

    %% defining and removing IED outliers.
    % removing outliers based on speed.
	ampOutliers = range(IEDwf') > ampThresholds(2) | range(IEDwf') < ampThresholds(1);
	allOutliers = speedOutliers | ampOutliers;

    % removing outliers based on speed.
 %   IEDdateTimeObject(allOutliers) = [];
    IEDspeed(allOutliers) = [];
    IEDdirection(allOutliers) = [];
    sigIEDs(allOutliers) = [];
    
    %% summary circ stats ||DIRECTION||
    circStats(pt) = circ_stats(IEDdirection(sigIEDs)');
%    [HR_p(pt), HR_t(pt)] = hrtest(IEDdirection(sigIEDs)');
%    [R_p(pt), R_z(pt)] = circ_otest(IEDdirection(sigIEDs)');
    k_concParam(pt) = circ_kappa(IEDdirection(sigIEDs)');
     
    % dealing with seizure data
	seizureExpansionFiles = dir(sprintf('/storage/selected_data/elliotWorking/Data/IEDs/seizureData/%s*IWbetas.mat',upper(patientID)));
	% parsing seizure times
    for sz = 1:length(seizureExpansionFiles)
        % loading ictal wavefront info
        load(fullfile(seizureExpansionFiles(sz).folder,seizureExpansionFiles(sz).name))
        
        saveDir = '/storage/selected_data/elliotWorking/Data/IEDs/seizureData';
        load(fullfile(saveDir,sprintf('%s-%d_IDdata.mat',upper(patientID),sz)))
        
        % ictal wavefront directions
		if pt=='c7'
			IWdirection = angle(IWdata.V(1) + IWdata.V(2)*sqrt(-1))-deg2rad(50);
		else
			IWdirection = angle(IWdata.V(1) + IWdata.V(2)*sqrt(-1));
		end
        % figuring out ID and IED distributions realtive to the direction of the IW.
        % INTERICTAL
        [IEDmtestP,IEDmu,IEDul,IEDll] = circ_mtest(unwrap(IEDdirection(sigIEDs)'-IWdirection),0)
        IEDmedtestP = circ_medtest(unwrap(IEDdirection(sigIEDs)'-IWdirection),0)
        
        % ICTAL
		if ~isfield(IDdata,'SecGenIdx')
			IDdata.SecGenIdx = ones(1,length(IDdata.sigModel));
		end
        [IDmtestP,IDmu,IDul,IDll] = circ_mtest(IDdata.burstDirections(IDdata.sigModel  & IDdata.SecGenIdx)'-IWdirection,0)
        IDmedtestP = circ_medtest(unwrap(angle(IDdata.burstDirections(IDdata.sigModel  & IDdata.SecGenIdx))'-IWdirection),0)
        
        
        %% THIS IS WHERE (a panel of) FIGURE 2 STARTS %%
        figure(2)
        
        % plotting direction histograms
        subplot(3,2,1)
        polarhistogram(angle(IDdata.burstDirections(IDdata.sigModel & IDdata.SecGenIdx))-IWdirection,18,'DisplayStyle','Stairs','Normalization','pdf')
        title('ICTAL discharges')
        
        subplot(3,2,2)
        polarhistogram(IEDdirection(sigIEDs)-IWdirection,18,'DisplayStyle','Stairs','Normalization','pdf')
        title('INTERICTAL discharges')
        
        subplot(3,2,3)
        ch = compass(IWdata.V(1) + sqrt(-1)*IWdata.V(2));
        title(sprintf('ictal wavefront direction for %s sz %d',patientID,sz))
        
        subplot(3,2,5)
        hold on
        if IEDmtestP
            text(0,5,sprintf('significant mean test for IEDs. mean:%.2f upper/lower:[%.2f %.2f]',IEDmu,IEDll,IEDul))
        else
            text(0,5,sprintf('INSIGNIFICANT mean test for IEDs. mean:%.2f upper/lower:[%.2f %.2f]',IEDmu,IEDll,IEDul))
        end
        if IDmtestP
            text(0,4,sprintf('significant mean test for IDs. mean:%.2f upper/lower:[%.2f %.2f]',IDmu,IDll,IDul))
        else
            text(0,4,sprintf('INSIGNIFICANT mean test for IDs. mean:%.2f upper/lower:[%.2f %.2f]',IDmu,IDll,IDul))
        end
        text(0,3,sprintf('ictal wavefront direction: %.2f',IWdirection))
        text(0,2,sprintf('INTERICTAL median test p-value:%.2f',IEDmedtestP))
        text(0,1,sprintf('ICTAL discharge median test p-value:%.2f',IDmedtestP))
        ylim([0 8])
        axis off
        hold off
        
        subplot(3,2,6)
        hold on
        line([1 1],[IEDll IEDul],'color',[0.5 0.5 0.5])
        scatter(1,IEDmu,30,[0.5 0.5 0.5],'o','filled')
        line([1 1],[IDll IDul],'color',[0 0 0])
        scatter(1,IDmu,30,[0 0 0],'o','filled')
        hold off
        xlim([0 2])
        ylim([-pi pi])
        
        
        % saving.
        halfMaximize(2,'left')
        print(2,sprintf('~/%s_%d_IEDdistributions_wIWdirections',patientID,sz),'-dpdf','-fillpage')
        
        close (2)
  

    end
    
end % looping over patients.aaaa
































