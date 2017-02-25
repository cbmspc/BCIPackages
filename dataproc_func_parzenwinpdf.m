

function [f] = dataproc_func_parzenwinpdf(test, train, littlesigma)


if size(test,2) ~= size(train,2)
    error('# columns of test and of train must be equal');
end

p = size(train,2);
x = test;
Ntrain = size(train,1);
f = 0;

% if strcmp(wfac,'k')
%     w = kde2d(train);
%     Sigma = diag(w.^2);
% elseif strcmp(wfac,'r')
%     w = dataproc_func_parzenwinwidth(train);
%     Sigma = w*w*eye(p);
% end

if ~exist('littlesigma','var') || isempty(littlesigma) || (isnumeric(littlesigma) && norm(littlesigma) == 0) % unspecified or zero
    if p == 1
        w = dataproc_func_kde1d(train);
    elseif p == 2
        w = dataproc_func_kde2d(train);
    else
        w = dataproc_func_parzenwinwidth(train);
        w = w*ones(p,1);
    end
    Sigma = diag(w.^2);
elseif size(littlesigma,1) == p && size(littlesigma,2) == p % valid covariance matrix
    Sigma = littlesigma;
elseif isvector(littlesigma) && length(littlesigma) == p % valid vector of sigmas
    Sigma = diag(littlesigma.^2);
elseif isscalar(littlesigma) % valid scalar bandwidth factor
    Sigma = (littlesigma.^2)*eye(p);
end

invSigma = inv(Sigma);
denominator = sqrt( (2*pi)^p * det(Sigma) );
for j = 1:Ntrain
    xdiff = x - ones(size(x,1),1) * train(j,:);
    f = f + exp(-1/2 * sum((xdiff * invSigma) .* xdiff, 2)) / denominator; %#ok<MINV>
end



% if length(who('w')) & length(w)
%     Sigma = w*w*eye(p);
%     invSigma = inv(Sigma);
%     denominator = sqrt( (2*pi)^p * det(Sigma) );
%     for j = 1:Ntrain
%         xdiff = x - ones(size(x,1),1) * train(j,:);
%         f = f + exp(-1/2 * sum((xdiff * invSigma) .* xdiff, 2)) / denominator;
%     end
% else
%     for j = 1:Ntrain
%         Sigma = cov(train(j,:));
%         invSigma = inv(Sigma);
%         denominator = sqrt( (2*pi)^p * det(Sigma) );
%         xdiff = x - ones(size(x,1),1) * train(j,:);
%         f = f + exp(-1/2 * sum((xdiff * invSigma) .* xdiff, 2)) / denominator;
%     end
% end
f = f ./ Ntrain;

f(f<=eps(0)) = eps(0);

