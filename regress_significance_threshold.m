% Calculates the critical R-square value based on experimental data
% When a R-square is at or above this value, it is significant given the
% experimental data (X = indep var, Y = dep var) and P-value thres P_Thres
%
% If Rsq_Observed is specified and P_Thres is blank, Rsq_Observed is
% compared and p-value is returned
%
% Method: 0 = Naive mode
%         1 = Match experimental data covariance
%         2 = Match experimental data auto-regression

function [Rsq_Thres, Rsq_List] = regress_significance_threshold (X, Y, MCmax, P_Thres, Method, Rsq_Observed)
if sum(abs(X(:,1) - 1)) == 0
    % X is already affine. Delete the first column and create own affine
    X = X(:,2:end);
    Method = abs(Method);
elseif sum(abs(X(:,end) - 1)) == 0
    % X is already affine. Delete the last column and create own affine
    X = X(:,1:end-1);
    Method = abs(Method);
end
if ~exist('MCmax','var') || isempty(MCmax)
    MCmax = 10000;
end
if ~exist('P_Thres','var')
    P_Thres = [];
end
if ~exist('Method','var') || isempty(Method)
    Method = 2; % Use the most extensive method by default
end
if ~exist('Rsq_Observed','var')
    Rsq_Observed = [];
end

if isempty(P_Thres) && isempty(Rsq_Observed)
    
end

affine = 1;
% if Method < 0
%     Method = -Method;
%     affine = 0;
% end
Method = round(Method);

f = ~any(~isfinite([X Y]),2);
X = X(f,:);
Y = Y(f,:);

[N, D] = size(X);
if Method == 2
    [Ad, Qd] = system_identify(X.');
elseif Method == 1
    try
        Xcc = chol(cov(X));
    catch
        % More dimensions than samples. Fall back to Method 2
        Method = 2;
        [Ad, Qd] = system_identify(X.');
    end
end
Rsq_List = nan(MCmax,1);

if Method == 0
    for i = 1:MCmax
        if affine
            Xnew = [ones(N,1) rand(N,D)];
        else
            Xnew = rand(N,D);
        end
        B = (Xnew'*Xnew)\(Xnew'*Y);
        Y_hat = Xnew*B;
        Rsq_List(i) = 1 - sum((Y - Y_hat).^2) / sum((Y - mean(Y)).^2);
    end
elseif Method == 1
    for i = 1:MCmax
        if affine
            Xnew = [ones(N,1) randn(N,D)*Xcc];
        else
            Xnew = randn(N,D)*Xcc;
        end
        B = (Xnew'*Xnew)\(Xnew'*Y);
        Y_hat = Xnew*B;
        Rsq_List(i) = 1 - sum((Y - Y_hat).^2) / sum((Y - mean(Y)).^2);
    end
elseif Method == 2
    for i = 1:MCmax
        X_sim = zeros(D,N);
        for k = 2:N
            X_sim(:,k) = Ad*X_sim(:,k-1) + Qd*randn(D,1);
        end
        if affine
            Xnew = [ones(N,1) X_sim.'];
        else
            Xnew = X_sim.';
        end
        B = (Xnew'*Xnew)\(Xnew'*Y);
        Y_hat = Xnew*B;
        Rsq_List(i) = 1 - sum((Y - Y_hat).^2) / sum((Y - mean(Y)).^2);
    end
end

Rsq_List = sort(Rsq_List);

if ~isempty(P_Thres)
    for i = length(P_Thres):-1:1
        Rsq_Thres(i) = Rsq_List(end-floor(P_Thres(i)*MCmax)+1);
    end
elseif ~isempty(Rsq_Observed)
    for i = length(Rsq_Observed):-1:1
        Rsq_Thres(i) = nnz(Rsq_List >= Rsq_Observed(i)) / length(Rsq_List);
    end
end

return


function [Ad, Qd] = system_identify (x)

x = detrend(x.').';

dim_x = size(x,1);

M = size(x,2);

a = zeros(dim_x);
for k = 2:M
    a = a + (x(:,k) * x(:,k-1)');
end

b = zeros(dim_x);
for k = 2:M
    b = b + (x(:,k-1) * x(:,k-1)');
end

Ad = a / b;


a = zeros(dim_x);
for k = 2:M
    a = a + x(:,k)*x(:,k)';
end

b = zeros(dim_x);
for k = 2:M
    b = b + x(:,k-1) * x(:,k)';
end

Qd = 1/(M-1)*(a - Ad*b);


