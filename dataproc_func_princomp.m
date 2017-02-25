function [Un, Eval] = dataproc_func_princomp (data)
% function [Un, Eval] = princomp (data)
% Principal Component Analysis of data (n observations by p dimensions).
% Eval = Eigenvalues
% Un = Eigenvectors


[n p] = size(data);
dmdata = data - ones(n,1)*mean(data,1);


% % SVD Method
% [U S V] = svd(dmdata,'econ');
% s = diag(S);
% Eval = s(1:end-(n<=p)).^2/(n-1);
% Un = V(:,1:length(Eval));


% EIG Method
if p >= n  % Small sample size
    [V D] = eig(dmdata*dmdata');
    % % LAPACK Approximation
    % OPTS.issym = true;
    % OPTS.disp = 0;
    %[V D] = eigs(dmdata*dmdata',n-1,'LA',OPTS);
    [d,I] = sort(diag(D));
    I = flipud(I);
    d = d(I(1:n-1));
    Eval = d/(n-1);
    Un = V(:,I(1:n-1));
    Un = dmdata' * Un * diag(d.^(-0.5));
else  % Small dimension
    [Un D] = eig(dmdata'*dmdata);
    Eval = diag(D)/(n-1);
    [Eval,I] = sort(Eval);
    Eval = flipud(Eval);
    Un = Un(:,flipud(I));
    Eval(Eval<eps) = 0;
end

