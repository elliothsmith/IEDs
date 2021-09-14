function [p_val,k_stat,k_perm] = evaluateNormality(angles, nPerms, permSize, fName)

% this function will fit a VonMises distribution to the input angles
% and return a measure of how well the epirical distribution fits, 
% based on a Kuiper test that evaluates signifcance using [nPerms] permutations
% of angles in groups of [permSize].
% 
% adding a fourth input argument will save an anaylsis figure with that file name
%

% fitting the VonMises distribtution
[thetahat kappa] = circ_vmpar(angles);

% now sampling from that distribution
nAngles = length(angles);
alpha = circ_vmrnd(thetahat, kappa, nAngles)

% running the kuiper test between the original and ideal distributions
fprintf('\nrunning permutation kuiper test. This may take a moment.')
[p_val,k_stat,k_perm] = permutation_kuipertest(angles, alpha, nPerms, permSize,90,0);


if nargin>3 
	figure(1)

	subplot(3,2,1)
	polarhistogram(angles,18,'DisplayStyle','stairs')
	title('histogram of input angles')

	subplot(3,2,2)
	polarhistogram(alpha,18,'DisplayStyle','stairs')
	title('histogram of generated angles')

	subplot(3,2,[3 4])
	H = histfit(k_perm);
	xlabel('permutation statistics')
	title(sprintf('true data test stat = %.2f',k_stat))

	saveas(1,fName)
	close(1)
	
end




