% Whitening transformation

function [Aw, W] = whitening (x)
[V, D] = eig(cov(x,1));
W = (sqrt(D) \ (V'))';
Aw = x * W;
Aw = Aw - ones(size(Aw,1),1) * mean(Aw,1);

