function f = logmat(M)
% Logarithm of Matrix

%    f = logm(M);

    [V D] = eig(M);
    Lambda = diag(D)+eps;
    f = V*diag(log(Lambda))*V'; 

    return;

