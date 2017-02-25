% Whitening transformation
% W is the transformation matrix that causes whitening
% X is the whitened version of x
% M is the mean of X
%
% To reverse the process, x = (X+ones(size(X,1),1)*M) * inv(W)

function [W X M] = dataproc_func_whiten (x)
[V D] = eig(cov(x,1));
W = (sqrt(D) \ (V'))';
X = x * W;
M = mean(X,1);
X = X - ones(size(X,1),1) * M;
