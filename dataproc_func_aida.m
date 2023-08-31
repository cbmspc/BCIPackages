
% Feature extraction based on Approximated Information Discriminant
% Analysis by Zoran Nenadic

% TrainData is (obs x dim) (Each row is an observation)


function Fmat = dataproc_func_aida(TrainData, TrainLabels, MfeatureSpace)

if any(isnan(TrainData(:)))
    error('TrainData contains NaN');
end

if any(isnan(TrainLabels(:)))
    error('TrainLabels contains NaN');
end


if size(TrainData,2) <= MfeatureSpace
    % There is no valid reason to do feature extraction if input dimension
    % is already equal to output dimension
    Fmat = eye(size(TrainData,2));
    return
end

%
% Aida = [];
% mu_ida = [];
% Alda = [];
% mu_lda = [];
% Anew = [];
% mu_new = [];
% Ftrain = [];
% CrossValidError_l = [];
% CrossValidError_q = [];
% NdataSpace = [];


BinData = TrainData .';
classes = unique(TrainLabels);
Nclass = length(classes);
Nparm = size(TrainData,2);

Memship = zeros(length(TrainLabels),1);

for i = 1:Nclass
    Memship(TrainLabels == classes(i)) = i-1;
end

% Tolx = 10^(-10);
% Tolf = 10^(-10);
% Niter = 5;

% Create an optimization options structure with optimization function
% FMINUNC that finds the minimum of a function of several variables.

% OPT = optimset('fminunc');
% OPT = optimset(OPT,'GradObj','on','Hessian','on','Display','off', ...
%     'TolFun',Tolf, 'TolX',Tolx,'MaxIter',2000);

% Here, unlike extract_ida_features.m, we do not set aside one trial/observation
% as validatory, and do not compute a number of sum(Ntr) transformation
% matrices for each observation excluded as validatory.
% We compute just one transformation matrix, and validate classifier later
% in on-line experiment.

%need to estimate class-conditional means and covariances

n = size(BinData,1);
m = MfeatureSpace;

M = zeros(n,1);
Sw = zeros(n,n);
Sb = zeros(n,n);
Si = zeros(Nparm,Nparm,Nclass);

p_i = zeros(1,Nclass);
Mi = zeros(Nparm,Nclass);

for j = 1:Nclass   % Left : Right :Up :Down
    lind = Memship == j-1;
    nind = nnz(lind);
    p_i(j) = nind/length(Memship);
    Mi(:,j) = mean(BinData(:,lind),2);
    Si(:,:,j) = cov(BinData(:,lind)');
    if nind > 1
        Si(:,:,j) = cov(BinData(:,lind)');
    else
        Si(:,:,j) = zeros(Nparm,Nparm);
    end


    %overall mean
    M = M + p_i(j) * Mi(:,j);

    %within class matrix
    Sw = Sw + p_i(j) * Si(:,:,j);

    %between class matrix
    Sb = Sb + p_i(j) * (Mi(:,j) * Mi(:,j)');
end

Sb = Sb - M * M';
S = Sw + Sb;

%pause

%display options for eigs
%OPTS.disp = 0;

%------------------- LDA (Linear Discriminant Analysis)

% if m < 2
%     [V D]= eigs(Sb,Sw,m,'LM',OPTS);
%     Alda = V/norm(V);
%     Alda = Alda';
%     mu_lda = -dataproc_func_get_negative_mu(Alda,S,Si,p_i);
% else
%     Alda = [];
%     mu_lda = 0;
% end

%------------------- NEW ALGORITHM (USE AS INITIAL CONDITION)

rootSw = dataproc_func_sqrtmat(Sw);   %transforms the variables so Sw = I;
W = inv(rootSw);      %whitening matrix so that Sw = I;

Znew = zeros(n,n);

Zi = zeros(Nparm,Nparm,Nclass);
for k = 1:Nclass
    Zi(:,:,k) = W * Si(:,:,k) * W;  %#ok<MINV> %transformed si's
    Znew = Znew - p_i(k) * dataproc_func_logmat(Zi(:,:,k));
end
Znew = dataproc_func_logmat(W*S*W) + Znew; %#ok<MINV>
[V, D] = eig(Znew);
[~, iD] = sort(diag(D));
V = fliplr(V(:,iD));
V = V(:,1:m);
Anew = V'*W; %#ok<MINV>
Anew = Anew/norm(Anew);
%mu_new = -dataproc_func_get_negative_mu(Anew,S,Si,p_i);


Fmat = Anew.';

