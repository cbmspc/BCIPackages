function [f, logf] = dataproc_func_normalpdf(test, train, Sigma)

m1 = size(test,2);
m2 = size(train,2);
if m1 ~= m2
    %warning('# columns of test (%i) and of train (%i) should be equal',m1,m2);
    m = min([m1 m2]);
    test = test(:,1:m);
    train = train(:,1:m);
end

if isempty(who('Sigma')) || isempty(Sigma)
    Sigma = cov(train);
end

xdiff = test - ones(size(test,1),1) * mean(train,1);
%f = exp(-1/2 * sum((xdiff * inv(Sigma)) .* xdiff, 2)) / sqrt( (2*pi)^size(train,2) * det(Sigma) );
f = exp(-1/2 * sum((xdiff/Sigma) .* xdiff, 2)) / sqrt( (2*pi)^size(train,2) * det(Sigma) );

sqrt_2pi_term = size(train,2)*log(2*pi)/2;
logf = -1/2 * sum((xdiff/Sigma) .* xdiff, 2) -log(sqrt(det(Sigma))) -sqrt_2pi_term;
