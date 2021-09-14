function Hx = XH(P,Q)
%   Hx = XH(P,Q) computes the cross entropy of two distributions. 
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
		Hx = - nansum( P .* log2( Q ) );

