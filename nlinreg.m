% Non-linear least-squares regression fit
% Example modelfun = @(beta,x)beta(1)+x*beta(2:end)
function [B, Covar, Resid, rmse, Y_hat, COD, ACOD, PVAL] = nlinreg (X, Y, modelfun, cv)

if ~exist('cv','var') || isempty(cv)
    cv = 0;
end

X = abs(X) + eps;

n = size(X,1);

B0 = rand(1+size(X,2),1);
opts = statset('nlinfit');
opts.RobustWgtFun = 'bisquare';

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
    f = ~any([isnan(X) isnan(Y)],2);
    COD = nan(1,size(Y,2));
    for v = 1:size(Y,2)
        c = corrcoef(Y_hat(f,v), Y(f,v));
        COD(1,v) = c(1,2)^2;
    end
else
    B = nan(size(B0,1),size(Y,2));
    Resid = nan(size(Y));
    Y_hat = nan(size(Y));
    rmse = nan(1,size(Y,2));
    parfor v = 1:size(Y,2)
        [B(:,v), Resid(:,v), ~, ~, MSE] = nlinfit(X, Y(:,v), modelfun, B0, opts);
        Y_hat(:,v) = modelfun(B(:,v), X);
        rmse(v) = sqrt(MSE);
    end
    SS_res = nansum(Resid.^2); %#ok<*NANSUM>
    SS_tot = nansum((Y-ones(size(Y,1),1)*mean(Y,1)).^2);
    COD = 1 - SS_res./SS_tot;
end
Covar = Resid'*Resid/size(Y,1);
ACOD = 1 - (1 - COD) * (n-1)/(n-p-1);

PVAL = nan(1,size(Y,2));
for v = 1:size(Y,2)
    [~,pv] = corrcoef(Y_hat(f,v), Y(f,v));
    PVAL(:,v) = pv(1,2).^2;
end

return


