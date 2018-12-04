% Finds the best lag resulting in highest positive covariance
%
% Calculate the normalized cross covariance between sig1 and sig2, and
% returns the coefficient (c), Student's T-test p-value (p), and the lag
% with the highest coefficient. To get goodness-of-fit, square c.
% 
% shiftlim limits the allowable lead/lag shift.
%
% Note: If dimensions are different, the longer signal will be truncated at
% the end
%
% SWliteral: Calculates the Pearson product-moment correlation coefficient
% literally without using xcov.

function [c, p, lag, rmse1, rmse2, lcpair] = covcoef (sig1, sig2, LagLimit, SWliteral)
if ~exist('LagLimit','var') || isempty(LagLimit)
    LagLimit = 0;
end
if ~exist('SWliteral','var') || isempty(SWliteral)
    SWliteral = 0;
end

if size(sig1,1) == 1 && size(sig1,2) > 1
    sig1 = sig1.';
end
if size(sig2,1) == 1 && size(sig2,2) > 1
    sig2 = sig2.';
end

d1 = size(sig1,2);
d2 = size(sig2,2);
l1 = size(sig1,1);
l2 = size(sig2,1);

if exist('LagLimit', 'var') && ~isempty(LagLimit)
    WarnLen = max(abs(LagLimit)) * 2;
else
    WarnLen = 10;
end
if l1 <= 0 || l2 <= 0
    error('An input signal is empty (zero length).');
end
if l1 < WarnLen || l2 < WarnLen
    warning('The input signals for covcoef are really short. Warning threshold = %i. Sig1_Len = %i. Sig2_Len = %i.', WarnLen, l1, l2);
end
c = nan(d1, d2);
p = c;
lag = c;
rmse1 = c;
rmse2 = c;
lcpair = cell(d1, d2);
for i = 1:size(sig1,2)
    for j = 1:size(sig2,2)
        [c(i,j), p(i,j), lag(i,j), rmse1(i,j), rmse2(i,j), lcpair{i,j}] = covcoef_func1 (sig1(:,i), sig2(:,j), LagLimit, SWliteral);
    end
end
if d1 == 1 && d2 == 1
    lcpair = lcpair{1};
end
return


function [c, p, lag, rmse1, rmse2, lcpair] = covcoef_func1 (sig1, sig2, LagLimit, SWliteral)

if ~exist('SWliteral','var') || isempty(SWliteral)
    SWliteral = 0;
end

sig1 = sig1(:);
sig2 = sig2(:);

L1 = length(sig1);
L2 = length(sig2);
if L2 > L1
    sig2 = sig2(1:L1);
elseif L1 > L2
    sig1 = sig1(1:L2);
end

nini = ~isnan(sig1) & ~isnan(sig2);
sig1 = sig1(nini);
sig2 = sig2(nini);

if size(sig1,1) == 0 || size(sig2,1) == 0
    c = NaN;
    p = NaN;
    lag = NaN;
    rmse1 = NaN;
    rmse2 = NaN;
    lcpair = NaN;
    return
end

sig1 = sig1-mean(sig1);
sig2 = sig2-mean(sig2);

if exist('LagLimit','var') && ~isempty(LagLimit)
    if length(LagLimit) == 1
        LagLimit = abs(LagLimit);
        LagLimit = [-1 1]*LagLimit;
    end
    maxlag = max(abs(LagLimit));
else
    maxlag = [];
end


if SWliteral
    % Do cross-covariance "manually"
    Lags = LagLimit(1):LagLimit(2);
    c = zeros(1,length(Lags));
    lags = c;
    parfor i = 1:length(Lags)
        lag = round(Lags(i));
        
        if lag > 0
            s1 = sig1(1+lag:end);
            s2 = sig2(1:end-lag);
        elseif lag < 0
            s1 = sig1(1:end+lag);
            s2 = sig2(1-lag:end);
        else
            s1 = sig1;
            s2 = sig2;
        end
        
        % Add a gradient to get rid of potential NaNs
        s1 = s1 + getepsgrad(s1);
        s2 = s2 + getepsgrad(s2);
        
        %cv = cov(s1, s2);
        %c(i) = cv(1,2) / (std(s1)*std(s2));
        
        % 20120912: If the length is exactly one, corrcoef does not work.
        if length(s1) > 1 && length(s2) > 1
            r = corrcoef(s1,s2);
        else
            r = eye(2);
        end
        
        c(i) = r(1,2);
        lags(i) = lag;
    end
    
else
    
    if ~isempty(maxlag)
        [c lags] = xcov(sig1,sig2,maxlag,'coeff');
    else
        [c lags] = xcov(sig1,sig2,'coeff');
    end
    
end

c(lags<min(LagLimit)) = -inf;
c(lags>max(LagLimit)) = -inf;

lcpair = [lags(:), c(:)];

[c i] = max(c);
lag = round(lags(i));

if lag > 0
    s1 = sig1(1+lag:end);
    s2 = sig2(1:end-lag);
elseif lag < 0
    s1 = sig1(1:end+lag);
    s2 = sig2(1-lag:end);
else
    s1 = sig1;
    s2 = sig2;
end
rmse1 = mean((s1 - s2).^2)^(0.5);
s2b = equalize_traces(s1, s2);
rmse2 = mean((s1 - s2b).^2)^(0.5);

if SWliteral
    s1 = s1 + getepsgrad(s1);
    s2 = s2 + getepsgrad(s2);
end
[~, p] = corrcoef(s1, s2);
p = p(1,2);


% % Use affine linear regression to calculate p-value.
% Only works at zero-lag
% [~,~,~,~,S] = regress(sig1, [ones(length(sig2),1) sig2]);
% p = S(3);



function y = getepsgrad (x)
if ~isempty(x)
    s = size(x,1);
    y = (1:s).';
    y = y - mean(y);
    y = y * eps(max(y));
else
    y = 0;
end