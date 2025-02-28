% Whitening transformation
% x should be tall, where each row is a time sample
% OPTIONAL: If x_baseline is passed in as argument, use it instead of x to
% generate the W matrix. x_baseline should also be tall.
%
% OPTIONAL: Set Flag_sym = true to create a symmetric W matrix.
%   There are more than one ways to whiten. This is another way where W is
%   symmetric.
%
% Output:
% W is the transformation matrix that causes whitening. 
%    W post-multiplies x.
% X is the whitened version of x
% M is the mean of X
%
% To reverse the process, x = X * inv(W) + M
%

function [W, X, M] = dataproc_func_whiten (x, x_baseline, Flag_sym)
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

if exist('Flag_sym','var') && isscalar(Flag_sym) && Flag_sym
    W = W * (V');
end

X = x * W;
