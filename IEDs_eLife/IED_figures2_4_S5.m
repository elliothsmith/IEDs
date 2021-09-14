Fs = 3e4;
% looping over patients. 
for pt = 1:7

	clearvars -except Fs pt HR_p HR_t circStats R_p R_z k_concParam
	
	switch pt
		case 1
			patientID = '';
			startSz = 1;
	end

	% counting IEDs starting from 0
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

			% IED waveforms. 
			IEDwf(nIEDs + id,:) = thisIED(id).medianWF;

		end % looping over IEDs in each file.
		nIEDs = nIEDs + nIEDsFl;
	end % looping over files.

	% removing outliers based on speed. 
	speedOutliers = outliers(IEDspeed);
	ampOutliers = range(IEDwf') > 2*IQR*range(IEDwf');
	allOutliers = speedOutliers | ampOutliers;

% 	IEDdateTimeObject(speedOutliers) = [];
	IEDspeed(allOutliers) = [];
	IEDdirection(allOutliers) = [];
	sigIEDs(allOutliers) = [];

	%% summary circ stats
	circStats(pt) = circ_stats(IEDdirection(sigIEDs)');
	[HR_p(pt), HR_t(pt)] = hrtest(IEDdirection(sigIEDs)')
	[R_p(pt), R_z(pt)] = circ_otest(IEDdirection(sigIEDs)')
	k_concParam(pt) = circ_kappa(IEDdirection(sigIEDs)')

	% loading seizure data
	seizureExpansionFiles = dir(sprintf('/storage/selected_data/elliotWorking/Data/IEDs/seizureData/%s*IWbetas.mat',upper(patientID)));

  	% parsing seizure times
  	for sz = startSz:length(seizureExpansionFiles)
		% loading ictal wavefront info. 
		load(fullfile(seizureExpansionFiles(sz).folder,seizureExpansionFiles(sz).name))

		% loading ictal discharge data. 
		saveDir = '/storage/selected_data/elliotWorking/Data/IEDs/seizureData';
		load(fullfile(saveDir,sprintf('%s-%d_IDdata.mat',upper(patientID),sz)))

		% calculating ictal wavefront directions.
		IWdirection = unwrap(angle(IWdata.V(1) + IWdata.V(2)*sqrt(-1)));

		% figuring out ID and IED distributions realtive to the direction of the IW.
		% INTERICTAL
		[IEDmtestP IEDmu IEDul IEDll] = circ_mtest(unwrap(IEDdirection(sigIEDs)'),0)
		IEDmedtestP = circ_medtest(unwrap(IEDdirection(sigIEDs)'),0)
		% ICTAL
		[IDmtestP IDmu IDul IDll] = circ_mtest(unwrap(angle(IDdata.burstDirections(IDdata.sigModel))'),0)
		IDmedtestP = circ_medtest(unwrap(angle(IDdata.burstDirections(IDdata.sigModel))'),0)

		%% THIS IS WHERE (a panel of) FIGURE 2 STARTS %%
		figure(2)

		% plotting direction histograms
		subplot(3,2,1)
		polarhistogram(angle(IDdata.burstDirections(IDdata.sigModel)),18,'DisplayStyle','Stairs','Normalization','pdf')
		title('ICTAL discharges')

		subplot(3,2,3)
		polarhistogram(IEDdirection(sigIEDs),18,'DisplayStyle','Stairs')
		title('INTERICTAL discharges')
			
		subplot(3,2,5)
		compass(IWdata.V(1) + sqrt(-1)*IWdata.V(2));
		title(sprintf('ictal wavefront and ictal discharge directions for %s sz %d',patientID,sz))
		
		subplot(3,2,2)
		hold on
		line([0.8 1.2],[IWdirection IWdirection],'color',[0.1 0.1 0.5])
		line([1 1],[IEDll IEDul],'color',[0.5 0.5 0.5]) 
		scatter(1,IEDmu,30,[0.5 0.5 0.5],'o','filled')
		line([1 1],[IDll IDul],'color',[0 0 0]) 
		scatter(1,IDmu,30,[0 0 0],'o','filled')
		hold off
		xlim([0 2])
		ylim([-pi pi])	
		title('IW, ID and IED directions')

		subplot(3,2,4)
		hold on
		line([1 1],abs([IEDll IEDul]-repmat(IWdirection,1,2)),'color',[0.5 0.5 0.5]) 
		scatter(1,abs(IEDmu-IWdirection),30,[0.5 0.5 0.5],'o','filled')
		line([1 1],abs([IDll IDul]-repmat(IWdirection,1,2)),'color',[0 0 0]) 
		scatter(1,abs(IDmu-IWdirection),30,[0 0 0],'o','filled')
		hold off
		xlim([0 2])
		ylim([0 pi])	
		title('ID and IED directions i.r.t. IW direction.')

		subplot(3,2,6)
		hold on
		line([1 1],abs([IEDll IEDul]-repmat(IDmu,1,2)),'color',[0.5 0.5 0.5]) 
		scatter(1,abs(IEDmu-IDmu),30,[0.5 0.5 0.5],'o','filled')
		hold off
		xlim([0 2])
		ylim([0 pi])	
		title('IED directions i.r.t. IDs" mean direction.')

		% saving. 
		halfMaximize(2,'left')
		print(2,sprintf('~/%s_%d_distSummaryPlot',patientID,sz),'-dpdf','-fillpage')
		close (2)

 	end

end % looping over patients. 





