% Multivariable and multivariate least-squares regression to find optimal
% B such that the 2-norm of (Y - f(B,X)) is minimized. Both linear and
% non-linear versions are supported.
%
% SYNTAX:
% [B, Covar, Resid, rmse, Y_hat, COD, ACOD, N, pvalue] = mvlsq (X, Y, affine, fast, cv)
%
% Required Inputs:
% X (n by p): Independent variable with n samples and p dimensions. 
%             NaN samples are ignored.
% Y (n by d): Dependent variable with n samples and d dimensions. 
%             NaN samples are ignored.
%
% Optional Inputs:
% affine: 0 = Do not add a column of ones. (default)
%         1 = Add a column of ones as the coefficient term.
%         function_handle = Use this function for non-linear regression.
%         Example: affine = @(beta,x)beta(1)+x*beta(2:end)
% 
% fast: 0 = Enable robust fitting. (default)
%       1 = Disable robust fitting.
%
% cv: 0 = Return Covar, Resid, rmse, Y_hat, COD, and ACOD using default
%         least-squares method. (default)
%     1 = Use leave-one-out cross-validation for calculations of Covar,
%         Resid, rmse, Y_hat, COD, and ACOD. This can mitigate overfitting
%         during model selection.
%         
% 
%
% Outputs:
% B (p by d): Parameters matrix
%
% COD (1 by d): Coefficient of determination (R^2 value)
%
% ACOD (1 by d): Adjusted COD (R bar squared)
%
% N (1 by d): Number of finite (non-nan) samples
%
% pvalue: P-value of COD
%
% EXAMPLES:
% [B, ~, ~, ~, Y_hat, Rsq] = mvlsq (X, Y, 1);
%
% modelfun = @(beta,x)beta(1)+x*beta(2:end);
% [B2, ~, ~, ~, Y2_hat, Rsq2] = mvlsq (X2, Y2, modelfun);
% 
%
% OVERLOADED FUNCTION SYNTAX:
% Y_hat = mvlsq (X, B, affine)
% 

function [B, Covar, Resid, rmse, Y_hat, COD, ACOD, N, pvalue] = mvlsq (X, Y, affine, fast, cv)

if ~exist('affine','var') || isempty(affine)
    affine = 0;
end

if ~exist('fast','var') || isempty(fast)
    fast = 0;
end

if ~exist('cv','var') || isempty(cv)
    cv = 0;
end

n = size(X,1);
d = size(X,2);

if d > 10*n
    fprintf('mvlsq warning: Your data has %i samples and %i dimensions. Press Ctrl+C to abort.', n, d);
    pause(5.0);
    fprintf('\n\n');
end

if n <= 1
    error('Cannot work with only 1 sample');
end


if isa(affine,'function_handle')
    if size(Y,1) ~= n && nargout == 1
        % Overloaded syntax mode
        B = affine(Y,X);
        return
    end
    
    [B, Covar, Resid, rmse, Y_hat, COD, ACOD, N] = nlinreg(X, Y, affine, fast, cv);
    return
end

if affine
    X = [ones(n,1) X];
end

p = size(X,2);

if size(Y,1) ~= n && size(Y,1) == p && nargout == 1
    % Overloaded syntax mode
    B = X*Y;
    return
end



for i = size(Y,2):-1:1
    f(:,i) = ~any([~isfinite(X) ~isfinite(Y(:,i))],2);
    N(:,i) = nnz(f(:,i));
end

try
    for i = size(Y,2):-1:1
        [B(:,i), rmse(:,i)] = calc_coef(X(f(:,i),:), Y(f(:,i),i), fast, f);
    end
catch
    % Error in regression. Return nan values.
    B = nan(size(X,2),size(Y,2));
    rmse = nan(1,size(Y,2));
end

try
    if cv
        Y_hat = nan(size(Y));
        for i = 1:size(X,1)
            b_train = calc_coef(X([1:i-1,i+1:end],:), Y([1:i-1,i+1:end],:), fast, f([1:i-1,i+1:end],:)); 
            Y_hat(i,:) = X(i,:)*b_train;
        end
    else
        Y_hat = X*B;
    end
catch
    Y_hat = nan(size(Y));
end
Resid = Y - Y_hat;
%Covar = Resid'*Resid/size(Y,1);
Covar = cov(Resid);
if fast || cv || any(isnan(rmse))
    rmse = diag(Covar).';
end

% if cv
%     COD = nan(1,size(Y,2));
%     for i = 1:size(Y,2)
%         if p - logical(affine) == 0
%             COD(:,i) = NaN;
%         else
%             c = corrcoef(Y_hat(f(:,i),i), Y(f(:,i),i));
%             COD(:,i) = c(1,2).^2;
%         end
%     end
% else
%     for i = size(Y,2):-1:1
%         SS_res = sum(Resid(f(:,i),i).^2);
%         SS_tot = sum((Y(f(:,i),i) - mean(Y(f(:,i),i))).^2);
%         COD(1,i) = 1 - SS_res./SS_tot;
%     end
% end

COD = nan(1,size(Y,2));
for i = 1:size(Y,2)
    if p - logical(affine) == 0
        COD(:,i) = NaN;
    else
        c = corrcoef(Y_hat(f(:,i),i), Y(f(:,i),i));
        if numel(c) == 4
            COD(:,i) = c(1,2).^2;
        else
            COD(:,i) = NaN;
        end
    end
end

ACOD = 1 - (1 - COD) * (n-1)/(n-p-1);

warning('off', 'MATLAB:nearlySingularMatrix');
punder = zeros(1,size(Y,2));
if nargout >= 9
    % SigThres is specified. Use better estimates
    MCmax = 10000;
    Met = 2;
elseif nargout == 0
    MCmax = 1000;
    Met = 2;
else
    MCmax = 0;
    Met = 2;
end
if MCmax > 0
    for i = size(Y,2):-1:1
        pvalue(i) = regress_significance_threshold(X, Y(:,i), MCmax, [], Met, COD(i));
        if pvalue(i) < 1/MCmax, punder(i) = 1/MCmax; end
    end
end

if nargout == 0
    pvdisp = pvalue;
    for i = 1:size(Y,2)
        le = '=';
        if punder(i) > 0
            le = '<';
            pvdisp(i) = punder(i);
        end
        ColName = '';
        if p > 1
            ColName = sprintf('Col#%i',i);
        end
        fprintf('%-8s COD = %.3f (p%s%.2g)\n', ColName, COD(i), le, pvdisp(i));
        fprintf('           n = %i\n', N(i));
        fprintf('           B = [%s]\n', num2str(B(:,i).',3));
    end
end

return


function [B, rmse] = calc_coef (X, Y, fast, f)
warning('off', 'stats:statrobustfit:IterationLimit');
B = zeros(size(X,2),size(Y,2));
rmse = nan(1,size(Y,2));
if fast
    B = (X'*X)\(X'*Y);
    if isnan(B(1))
        % Backup method if X is rank deficient
        for i = 1:size(Y,2)
            B(:,i) = pinv(X(f(:,i),:))*Y(f(:,i),i);
        end
    end
else
    for i = 1:size(Y,2)
        try
            [B(:,i), S] = robustfit(X,Y(:,i),[],[],'off');
            rmse(i) = S.s;
        catch
            % Use non-robust method as a backup
            B(:,i) = (X'*X)\(X'*Y(:,i));
            if isnan(B(1,i))
                % Backup method if X is rank deficient
                B(:,i) = pinv(X(f(:,i),:))*Y(f(:,i),i);
            end
            %B(:,i) = nan(size(X,2),size(Y,2));
        end
    end
end


% Non-linear least-squares regression fit
% Example modelfun = @(beta,x)beta(1)+x*beta(2:end)
function [B, Covar, Resid, rmse, Y_hat, COD, ACOD, N] = nlinreg (X, Y, modelfun, fast, cv)

if ~exist('fast','var') || isempty(fast)
    fast = 0;
end

if ~exist('cv','var') || isempty(cv)
    cv = 0;
end

[n,p] = size(X);
for i = size(Y,2):-1:1
    f(:,i) = ~any([~isfinite(X) ~isfinite(Y(:,i))],2);
    N(:,i) = nnz(f(:,i));
end

B0 = rand(1+size(X,2),1);
opts = statset('nlinfit');
if ~fast
    opts.RobustWgtFun = 'bisquare';
    opts.Robust = 'on';
end

if cv
    B = nan(size(B0,1),size(Y,2));
    for v = 1:size(Y,2)
        B(:,v) = nlinfit(X, Y(:,v), modelfun, B0, opts);
        Y_hat = nan(size(Y));
        for i = 1:n
            B_train = nlinfit(X([1:i-1,i+1:end],:), Y([1:i-1,i+1:end],v), modelfun, B0, opts); 
            Y_hat(i,v) = modelfun(B_train, X(i,:));
        end
    end
    Resid = Y - Y_hat;
else
    B = nan(size(B0,1),size(Y,2));
    Y_hat = nan(size(Y));
    rmse = nan(1,size(Y,2));
    for v = 1:size(Y,2)
        [B(:,v), ~, ~, ~, MSE] = nlinfit(X, Y(:,v), modelfun, B0, opts);
        Y_hat(:,v) = modelfun(B(:,v), X);
        rmse(v) = sqrt(MSE);
    end
    Resid = Y - Y_hat;
end
Covar = Resid'*Resid/size(Y,1);
if cv
    rmse = diag(Covar);
end

COD = nan(1,size(Y,2));
for v = 1:size(Y,2)
    c = corrcoef(Y_hat(f(:,v),v), Y(f(:,v),v));
    COD(1,v) = c(1,2)^2;
end

ACOD = 1 - (1 - COD) * (n-1)/(n-p-1);

return
