function sbest = dataproc_func_parzenwinwidth (X)

dd = squareform((pdist(X)));
dd(find(dd==0)) = NaN;
%dd = min(dd);
dd = reshape(dd,1,[]);

% Min = min(nndist);
% Max = max(nndist);

N = size(X,1);
d = size(X,2);

K = internalfunc_sturges(dd);

% Edges = logspace(log(Min)/log(10),log(Max)/log(10),K+1);
% Centers = diff(Edges);
% [Frior,Sig] = hist(nndist,Centers);

[Frior,s] = hist(dd,K);

Frior = log(Frior+eps);
Fikelihood = -N.*d./2.*log(2.*pi.*s.^2);
Fosterior = Frior + Fikelihood;

[Y,I] = max(Fosterior);

sbest = s(I);

function K = internalfunc_sturges(X)
X = unique(X);
N = length(X);
K = ceil(log(N)/log(2)+1);

