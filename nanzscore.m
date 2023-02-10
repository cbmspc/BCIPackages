function [z, mu, sigma] = nanzscore (x)
% https://www.mathworks.com/matlabcentral/answers/249566-zscore-a-matrix-with-nan
mu = mean(x,'omitnan');
sigma = std(x, 'omitnan');
z = bsxfun(@rdivide, bsxfun(@minus, x, mu), sigma);
