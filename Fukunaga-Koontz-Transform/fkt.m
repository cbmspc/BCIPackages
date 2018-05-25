function [Y1,Y2,T,Eval]=fkt(X1,X2)
% Input arguments
% X1:   N1 x L matrix of class 1 object data
% X2:   N2 x L matrix of class 2 object data
% Outputs
% Y1:   N1 x L transformed matrix of class 1 object data
% Y2:   N2 x L transformed matrix of class 2 object data
% T:    L x L transformation matrix
% Eval: a vector of the eigenvalues of R1

% calculate sample autocorrelation matrices
[N1,~]=size(X1);
[N2,~]=size(X2);
S1=(X1'*X1)/N1;
S2=(X2'*X2)/N2;

% calculate transform matrix
[Phi,D1]=eig(S1+S2);
W=Phi*sqrt(inv(D1));
[Psi,D]=eig(W'*S1*W);
T=W*Psi;
Eval=diag(D);

% transform
Y1=X1*T;
Y2=X2*T;