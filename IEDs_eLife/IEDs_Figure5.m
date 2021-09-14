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
for pt = 1:7
    clearvars -except Fs pt
    switch pt
        case 1
    end

    % initializing firing rates
	binWidth = 0.01;
	histEdges = 0:binWidth:1;

	% initializing number of IEDs
	nIEDs = 0;
	
    % loop over the files in the directory.
	dirList = dir(sprintf('/storage/selected_data/elliotWorking/Data/IEDs/allIEDs_%s/*.mat',patientID));
	for fl = 1:length(dirList)
        load(fullfile(dirList(fl).folder,dirList(fl).name))

 		% loop over files.
        nIEDsFl = length(thisIED);
		fprintf('\nprocessing preprocessed IED file %d of %d (contains %d IEDs)',fl,length(dirList),nIEDsFl)
        for id = 1:nIEDsFl
            % other features in the data we may want to visualaize.
            sigIEDs(nIEDs + id) = thisIED(id).wavePval<0.05;
            V = pinv(thisIED(id).waveBetas(1:2));
            IEDspeed(nIEDs + id) = sqrt(sum(V.^2));
            IEDdirection(nIEDs + id) = angle(V(1) + V(2)*sqrt(-1));
            clear V
            
	 		% getting morphology
	        IEDwf(nIEDs + id,:) = thisIED(id).medianWF;

	        % getting mean firing rates across channels for each IED.
			try
				N(nIEDs + id,:) = histcounts([thisIED(id).muaTimes{1:96}],histEdges);
			catch
				N(nIEDs + id,:) = histcounts([thisIED(id).muaTimes{:}],histEdges);
			end

        end % looping over IEDs in each file.
        nIEDs = nIEDs + nIEDsFl;
    end % looping over files.
    
    % removing outliers based on speed.
	allOutliers = speedOutliers | ampOutliers;

    % change of variables. Actually, I'll have to retain these indices
    % inorder to rotate back to the original distribution
    sigDirections = IEDdirection(~allOutliers)';
    sigSpeeds = IEDspeed(~allOutliers)';
    % sigDateTime = IEDdateTimeObject(~allOutliers)';	
 
	%% getting all the summary stats for overall IED direction distribution.
    circStats(pt) = circ_stats(sigDirections);
    k_concParam = circ_kappa(sigDirections);
    
	% evaluating VM fit for whole distribution
	nPerms = 1000;
	permSize = 60;
	[p_valAll,k_statAll,k_permsAll] = evaluateNormality(sigDirections, nPerms, permSize,sprintf('%s_overall_IEDdistFit.pdf',patientID));

    %% do gaussian mixture modeling here.
	threshold = [0.3 0.7];
	GMM = fitmvmdist(sigDirections,2);
	[clusters, nLogL, PostP] = cluster(GMM,sigDirections);
	idxBothSharedDiag = PostP(:,1)>=threshold(1) & PostP(:,1)<=threshold(2);
	numInBoth = sum(idxBothSharedDiag);

    % fitting VonMises Distribution for each cluster
    circStatsCluster1(pt) = circ_stats(sigDirections(clusters==1 | idxBothSharedDiag));
    circStatsCluster2(pt) = circ_stats(sigDirections(clusters==2 | idxBothSharedDiag));

    % evaluating VM fit for each gaussian mixture
    [p_val_1,k_stat_1,k_perm_1] = evaluateNormality(sigDirections(clusters==1 | idxBothSharedDiag), nPerms, permSize, sprintf('%s_cluster1_IEDdistFit.pdf',patientID));
    [p_val_2,k_stat_2,k_perm_2] = evaluateNormality(sigDirections(clusters==2 | idxBothSharedDiag), nPerms, permSize, sprintf('%s_cluster2_IEDdistFit.pdf',patientID));

	% saving cluster stats
	save(sprintf('~/%s_VMclusterStats.mat',patientID),'circStatsCluster1','p_val_1','k_stat_1','k_perm_1','circStatsCluster2','p_val_2','k_stat_2','k_perm_2','circStats','p_valAll','k_statAll','k_permsAll');

    plotDirectionInfo = true;
	plotClusterFeatures = true;
    if plotDirectionInfo
	    %% THIS IS WHERE FIGURE STARTS %%
	    figure(2)
	    subplot(3,2,1)
	    polarhistogram(sigDirections(clusters==1 | idxBothSharedDiag),18,'facecolor',rgb('springgreen'),'edgecolor','none')

		subplot(3,2,2)
        raincloud_plot(sigSpeeds(clusters==1 | idxBothSharedDiag), 'box_on', 1, 'color', rgb('springgreen'))
        title('cluster 1')

        subplot(3,2,3)
        polarhistogram(sigDirections(clusters==2 | idxBothSharedDiag),18,'facecolor',rgb('darksalmon'),'edgecolor','none')

		subplot(3,2,4)
        raincloud_plot(sigSpeeds(clusters==2 | idxBothSharedDiag), 'box_on', 1, 'color', rgb('darksalmon'))
        title('cluster 2')

		% testing for differences in speed... All significant [20200323]
        [clustP,~,clustStats] = ranksum(sigSpeeds(clusters==2 | idxBothSharedDiag),sigSpeeds(clusters==1 | idxBothSharedDiag));
        subplot(3,2,5)
        hold on
 
 		text(0,0.3,sprintf('IED distribution vs. VM distreibution with same parameters: stat = %.2f, p-value = %.2f',k_statAll,p_valAll))
        text(0,0.2,sprintf('difference in mean angles (rad) = %d',abs(circStatsCluster2(pt).mean - circStatsCluster1(pt).mean)));
        text(0,0.1,sprintf('mean speed for cluster 1: %.2f. mean speed for cluster 2: %.2f.',mean(sigSpeeds(clusters==1 | idxBothSharedDiag)),mean(sigSpeeds(clusters==2 | idxBothSharedDiag))));
        text(0,0,sprintf('median speed for cluster 1: %.2f. median speed for cluster 2: %.2f.',median(sigSpeeds(clusters==1 | idxBothSharedDiag)),median(sigSpeeds(clusters==2 | idxBothSharedDiag))));

		if clustP < 0.05
	        text(0,-0.1,sprintf('significant difference between clusters: U = %.2f, p = %.2f',clustStats.ranksum,clustP))
	    else
	        text(0,-0.1,sprintf('no significant difference in speed between clusters'))
        end

		text(0,-0.2,sprintf('cluster 1 (n = %d) distribution vs. VM distreibution with same parameters: stat = %.2f, p-value = %.2f',sum(clusters==1),k_stat_1,p_val_1))
        text(0,-0.3,sprintf('cluster 2 (n = %d) distribution vs. VM distreibution with same parameters: stat = %.2f, p-value = %.2f',sum(clusters==2),k_stat_2,p_val_2))
		text(0,-0.5,sprintf('%d discharge overlap between clusters',numInBoth))

		hold off
		axis off
	
		halfMaximize(2,'left')
		print(2,sprintf('~/%s_2IEDpotentialPopulations',patientID),'-dpdf','-fillpage')
		
		close (2)
 
	elseif plotClusterFeatures


		%% THIS IS WHERE FIGURE STARTS %%
		figure(2)

		subplot(3,2,1)
	    polarhistogram(sigDirections(clusters==1 | idxBothSharedDiag),18,'facecolor',rgb('springgreen'),'edgecolor','none')
		title(sprintf('cluster 1: N = %d',sum(clusters==1)))

		subplot(3,2,2)
		polarhistogram(sigDirections(clusters==2 | idxBothSharedDiag),18,'facecolor',rgb('springgreen'),'edgecolor','none')
		title(sprintf('cluster 2: N = %d',sum(clusters==2)))

		% morphology
		tSec = linspace(0,1,30001);
		subplot(3,2,3)
		hold on
		patch([tSec fliplr(tSec)], [median(IEDwf(clusters==1,:))-std(IEDwf(clusters==1,:))./sqrt(sum(clusters==1))  fliplr(median(IEDwf(clusters==1,:))+std(IEDwf(clusters==1,:))./sqrt(sum(clusters==1)))],rgb('seagreen'),'facealpha',0.5,'edgecolor','none')
		patch([tSec fliplr(tSec)], [median(IEDwf(clusters==2,:))-std(IEDwf(clusters==2,:))./sqrt(sum(clusters==2))  fliplr(median(IEDwf(clusters==2,:))+std(IEDwf(clusters==2,:))./sqrt(sum(clusters==1)))],rgb('steelblue'),'facealpha',0.5,'edgecolor','none')
		plot(tSec,median(IEDwf(clusters==2,:)),'color',rgb('steelblue'))
		plot(tSec,median(IEDwf(clusters==1,:)),'color',rgb('seagreen'))
%		text(0,20,sprintf('cluster 1: N = %d',sum(clusters==1)),'color',rgb('seagreen'))
%		text(0,40,sprintf('cluster 2: N = %d',sum(clusters==2)),'color',rgb('steelblue'))
		hold off
		axis tight square

		% firing rate
		tPSTH = histEdges(1:end-1)+binWidth/2;
		subplot(3,2,5)
		hold on
		patch([tPSTH fliplr(tPSTH)], [mean(N(clusters==1,:))-std(N(clusters==1,:))./sqrt(sum(clusters==1))  fliplr(mean(N(clusters==1,:))+std(N(clusters==1,:))./sqrt(sum(clusters==1)))],rgb('seagreen'),'facealpha',0.5,'edgecolor','none')
		patch([tPSTH fliplr(tPSTH)], [mean(N(clusters==2,:))-std(N(clusters==2,:))./sqrt(sum(clusters==2))  fliplr(mean(N(clusters==2,:))+std(N(clusters==2,:))./sqrt(sum(clusters==2)))],rgb('steelblue'),'facealpha',0.5,'edgecolor','none')
		plot(tPSTH,mean(N(clusters==2,:)),'color',rgb('steelblue'))
		plot(tPSTH,mean(N(clusters==1,:)),'color',rgb('seagreen'))
		hold off
		axis tight square

		% now looking at relationship between firing rate and IED amplitude for each cluster.
		LM1 = fitlm(range(IEDwf(clusters==1,:),2),mean(N(clusters==1,tPSTH>=0.45 & tPSTH<=0.55),2))
		subplot(3,2,4)
		scatter(range(IEDwf(clusters==1,:),2),mean(N(clusters==1,tPSTH>=0.45 & tPSTH<=0.55),2),4,rgb('steelblue'),'filled')
		title(sprintf('amp & Fr correlation: p = %.2f, R^2 = %.2f',LM1.Coefficients{2,4},LM1.Rsquared.Ordinary))

		LM2 = fitlm(range(IEDwf(clusters==2,:),2),mean(N(clusters==2,tPSTH>=0.45 & tPSTH<=0.55),2))
		subplot(3,2,6)
		scatter(range(IEDwf(clusters==2,:),2),mean(N(clusters==2,tPSTH>=0.45 & tPSTH<=0.55),2),4,rgb('seagreen'),'filled')
		title(sprintf('amp & Fr correlation: p = %.2f, R^2 = %.2f',LM2.Coefficients{2,4},LM2.Rsquared.Ordinary))

 
% 		subplot(3,2,3)
% 		raincloud_plot(range(IEDwf(clusters==1,:),2))
% 		title('cluster 1')
% 		subplot(3,2,5)
% 		raincloud_plot(range(IEDwf(clusters==2,:),2))
% 		title('cluster 2')
% 
% 		subplot(3,2,4)
% 		raincloud_plot(mean(N(clusters==1,tPSTH>=0.45 & tPSTH<=0.55),2))
% 		subplot(3,2,6)
% 		raincloud_plot(mean(N(clusters==2,tPSTH>=0.45 & tPSTH<=0.55),2))

		halfMaximize(2,'left')
		print(2,sprintf('~/%s_gaussianClusterFeatures',patientID),'-dpdf','-fillpage')

		close (2)

	end

end % looping over patients.



