function [p_val, k_stat, k_perm] = permutation_kuipertest(alpha1, alpha2, nPerms, permutationSize, res, vis_on)

% [pval, k, K] = permutation_kuipertest(alpha1, alpha2, nPerms, permuationSize, vis_on)
%
%   The Kuiper two-sample permutation test tests whether the two samples
%   differ significantly based on probabilities derived from a permutation
%   distribution. The difference can be in any property, such as mean
%   location and dispersion. It is a circular analogue of the
%   Kolmogorov-Smirnov test.
% 
%   H0: The two distributions are identical.
%   HA: The two distributions are different.
%
% Input: 
%   alpha1    fist sample (in radians)
%   alpha2    second sample (in radians)
%	permSize  The size of the samples drawn from each sampple without replacement.
%   res       resolution at which the cdf is evaluated
%   vis_on    display graph
%
% Output:
%   p_val       permutation p-value
%   k_stat      test statistic for real data
%   k_perm		test statistics for permutation data [length(k_perm) == nPerms]
%
% Modified from the Circular Statistics Toolbox for Matlab
% EHS20200325

% Update 2012
% By Marc J. Velasco and Philipp Berens, 2009
% velasco@ccs.fau.edu


if nargin < 4; 
	res = 90; 
elseif nargin < 5
    vis_on = 0;
end

n = length(alpha1(:));
m = length(alpha2(:));

% create cdfs of both samples
[phis1 cdf1 phiplot1 cdfplot1] = circ_samplecdf(datasample(alpha1,floor(permutationSize./2)), res);
[foo, cdf2 phiplot2 cdfplot2] = circ_samplecdf(datasample(alpha2,floor(permutationSize./2)), res); %#ok<ASGLU>

% maximal difference between sample cdfs
[dplus, gdpi] = max([0 cdf1-cdf2]);
[dminus, gdmi] = max([0 cdf2-cdf1]);

% calculate k-statistic
k_stat = sqrt(n + m) * (dplus + dminus);

% find p-value
% parfor_progress(nPerms);
k_perm = NaN(1,nPerms);
for ps = 1:nPerms
    perm_samp = datasample([alpha1(:); alpha2(:)],permutationSize,'Replace',false);
    
    % create cdfs of both samples
    [phis1 cdf1 phiplot1 cdfplot1] = circ_samplecdf(perm_samp(1:floor(permutationSize./2)), res);
    [foo, cdf2 phiplot2 cdfplot2] = circ_samplecdf(perm_samp(floor(permutationSize./2):end), res); %#ok<ASGLU>

    % maximal difference between sample cdfs
    [dplus, gdpi] = max([0 cdf1-cdf2]);
    [dminus, gdmi] = max([0 cdf2-cdf1]);

    % calculate k-statistic
    k_perm(ps) = sqrt(permutationSize) * (dplus + dminus);
%    parfor_progress
end
% parfor_progress(0);

% output values
p_val = sum(k_perm<k_stat)./nPerms;

% visualize
if vis_on
	% binning pemrutation statistics. 
	n_stat_bins = 100;
	[K_dist,perm_edges] = histcounts(k_perm,n_stat_bins);

	% plotting 
    histfit(gca,k_perm,n_stat_bins);
    xlabel('test statistic')
    ylabel('count')
    title(sprintf('true data test statistic = %d',k_stat))
end

end % eof
