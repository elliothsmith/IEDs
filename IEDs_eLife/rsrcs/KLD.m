function D = KLD(P,Q)
%      D = KLD(P,Q)
% 
%   Compute Kullback-Leibler divergence of probability distribution Q from probability distribution P.
%   P represents the "true" distribution of data, observations, or a theoretical distribution. 
%   whereas Q typically represents a theory, model, description, or approximation of P.
%   Both P and Q are thus distribution (obtained e.g. using 'getDensity.m') computed over the same range
%   (i.e. the same Bins -- if not, obviously it makes no sense to compare them).
% 
%   The Kullback-Leibler divergence is an non-symmetric measure (see below) of the difference between 
%   two probability distributions P and Q. Specifically, the Kullback–Leibler divergence of Q from P, 
%   is a measure of the information lost when Q is used to approximate P. The Kullback–Leibler divergence 
%   measures the expected number of extra bits (so intuitively it is non negative) required to code samples 
%   from P when using a code optimized for Q, rather than using the true code optimized for P. 
%   Although it is often intuited as a metric or distance, the Kullback–Leibler divergence is not a true 
%   metric — for example, it is not symmetric: the Kullback–Leibler divergence from P to Q is generally 
%   not the same as that from Q to P.
% 
%   NOTE:
%   The code treat P and Q as DISCRETE probability distributions
% 
%   Description inspired by Wikipedia page: https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence
%   ----------------------------------------------------------------------------------------------------------------------------------------------------------
% 
%
%   INPUT
% 
%           P    :   Reference Probability distribution (an <M x 1> vector, where M is the no. of bins used to compute the distribution)
%           Q    :   Probability distribution I want to compare with the reference one (an <M x 1> vector too)
%         
%   OUTPUT
% 
%           KLD  :   Kullback-Leibler divergence (a scalar)
%         
%         
%   --------------------------------------------------------------------------------------------------------------------------------------------------------
%   Ruggero G. Bettinardi, PhD student
% 
%   Computational Neuroscience Group,
%   Center for Brain & Cognition, Pompeu Fabra University
%   m: rug.bettinardi@gmail.com
%   --------------------------------------------------------------------------------------------------------------------------------------------------------
if abs(sum(P)-sum(Q)) > 1e-10
	    error('WARNING: Are you sure that P and Q are probability distributions ??? Apparently they do not sum up to 1 ... check : sum(P) & sum(Q)')
	end
	MP = numel(P);   % no. of bins of P
	MQ = numel(Q);   % no. of bins of Q
	if MP ~= MQ
		    error('WARNING: P and Q have different size!')
		end
		M   = numel(P);                   
		P   = reshape(P,[M,1]);           
		Q   = reshape(Q,[M,1]);
		D = nansum( P .* log2( P./Q ) );

