function f = sqrtmat(M)
% Square root of Matrix
%   f = sqrtm(M);
[V D] = eig(M);
Lambda = diag(D);
f = V*diag(sqrt(Lambda))*V';
return;
