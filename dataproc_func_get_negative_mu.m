function [mu, grad_mu, hess_mu] = get_negative_mu(T,S,Si,p_i)

%[mu, grad_mu, hess_mu] = get_negative_mu(T,S,Si,p_i)
%T - m x n transformation matrix
%S - ovarall covariances matrix
%Si - class conditional covariance matrices
%p_i - prior probabilities of classes


C = length(p_i);
[m n] = size(T);

%%GET THE FUNCTIONAL (NEGATIVE) MU
muo = 0;
for i = 1:C
  mu = muo + p_i(i) * log(det(T*Si(:,:,i)*T'));
  muo = mu;
end
mu = 0.5 * [log(det(T*S*T')) - muo];
%take the negative of mu
mu = -mu;

%%GET THE GRADIENT OF MU WITH RESPECT TO TRANSFORMATION MATRIX T
if nargout > 1
  grado = 0;
  for i = 1:C
    tmpmatr = T * Si(:,:,i) * T';
    %use ``pinv'' for (near) singular covariances 
    %grad_mu = grado + p_i(i) * pinv(tmpmatr) * T * Si(:,:,i);
    grad_mu = grado + p_i(i) * inv(tmpmatr) * T * Si(:,:,i);
    grado = grad_mu;
  end
  
  grad_mu = inv(T*S*T') * T * S - grado;
  grad_mu = -grad_mu;
end

%%GET THE HESSIAN OF MU WITH RESPECT TO TRANSFORMATION MATRIX T
if nargout > 2
    kk = 0;
    for i = 1:n
        for j = 1:m
            kk = kk + 1;
            b(kk) = kk;
            a(kk) = i + (j-1)*n;
        end
    end
    Th = sparse(a,b,ones(1,m*n));
    hesso = 0;
    for i = 1:C
        temp1= T * Si(:,:,i);
        %use ``pinv'' for near sinular covariances
        %temp2 = pinv(temp1 * T');
        temp2 = inv(temp1 * T');
        I = kron(eye(n),temp2)*kron(Si(:,:,i),eye(m))';
        II = -kron(temp1',eye(m))*kron(temp2',temp2) * ...
            [kron(temp1,eye(m)) + kron(eye(m),temp1)*Th];
        hess_mu = hesso + p_i(i) * [I + II];
        hesso = hess_mu;
    end
    temp1 = T * S;
    temp2 = inv(temp1 * T');
    I = kron(eye(n),temp2)*kron(S,eye(m))';        
    II = -kron(temp1',eye(m))*kron(temp2',temp2) * ... 
        [kron(temp1,eye(m)) + kron(eye(m),temp1)*Th];
    hess_mu = I + II - hesso;
    hess_mu = -hess_mu;
end