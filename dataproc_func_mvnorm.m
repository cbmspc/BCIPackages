% Multivariate normal distribution of training data X (obs x parm) at
% observation x. Rows of X are observations. Columns of X are parameters.
%

function f = mvnorm (X, x)

p = size(X,2);

Sigma = cov(X);

xdiff = x - ones(size(x,1),1) * mean(X,1);

f = exp(-1/2 * xdiff * inv(Sigma) * xdiff.') /sqrt((2*pi)^p * det(Sigma));

