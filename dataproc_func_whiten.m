% Whitening transformation
% W is the transformation matrix that causes whitening
% X is the whitened version of x
% M is the mean of X
%
% OPTIONAL: If x_baseline is passed in as argument, use it instead of x to
% generate the W matrix
%
% To reverse the process, x = X * inv(W) + M
%
% Outdated: To reverse the process, x = (X+ones(size(X,1),1)*M) * inv(W)

function [W, X, M] = dataproc_func_whiten (x, x_baseline)
M = mean(x,1);
x = x - ones(size(x,1),1)*M;

if exist('x_baseline','var')
    M_baseline = mean(x_baseline,1);
    x_baseline = x_baseline - ones(size(x_baseline,1),1)*M_baseline;
    [V, D] = eig(cov(x_baseline,1));
else
    [V, D] = eig(cov(x,1));
end

W = (sqrt(D) \ (V'))';
X = x * W;
