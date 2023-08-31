
% Feature extraction based on Information Discriminant Analysis by Zoran
% Nenadic, using Approximated IDA as initial condition.

% TrainData is (obs x dim) (Each row is an observation)

function Fmat = featureextraction(TrainData, TrainLabels, MfeatureSpace)
Aida = [];
mu_ida = [];
Alda = [];
mu_lda = [];
Anew = [];
mu_new = [];
Ftrain = [];
CrossValidError_l = [];
CrossValidError_q = [];
NdataSpace = [];

if any(isnan(TrainData(:)))
    error('TrainData contains NaN');
end

if any(isnan(TrainLabels(:)))
    error('TrainLabels contains NaN');
end


BinData = TrainData .';
classes = unique(TrainLabels);
Nclass = length(classes);
for i = 1:Nclass
    Memship(find(TrainLabels == classes(i))) = i-1;
end

Tolx = 10^(-10);
Tolf = 10^(-10);
Niter = 5;

% Create an optimization options structure with optimization function
% FMINUNC that finds the minimum of a function of several variables.
OPT = optimset('fminunc');
OPT = optimset(OPT,'GradObj','on','Hessian','on','Display','off', ...
    'TolFun',Tolf, 'TolX',Tolx,'MaxIter',2000);

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
for j = 1:Nclass   % Left : Right :Up :Down
    ind = find(Memship == j-1);    % Memship Left = 0, Right = 1
    p_i(j) = length(ind)/length(Memship);
    %p_i(j) = 0.5; % assume equal prior probability of L and R
    Mi(:,j) = mean(BinData(:,ind),2);
    Si(:,:,j) = cov(BinData(:,ind)');


    %overall mean
    M = M + p_i(j) * Mi(:,j);

    %within class matrix
    Sw = Sw + p_i(j) * Si(:,:,j);

    %between class matrix
    Sb = Sb + p_i(j) * [Mi(:,j) * Mi(:,j)'];
end

Sb = Sb - M * M';
S = Sw + Sb;

%pause

%display options for eigs
OPTS.disp = 0;

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

for k = 1:Nclass
    Zi(:,:,k) = W * Si(:,:,k) * W;  %transformed si's
    Znew = Znew - p_i(k) * dataproc_func_logmat(Zi(:,:,k));
end
Znew = dataproc_func_logmat(W*S*W) + Znew;
[V, D] = eig(Znew);
[D iD] = sort(diag(D));
V = fliplr(V(:,iD));
V = V(:,1:m);
Anew = V'*W;
Anew = Anew/norm(Anew);
mu_new = -dataproc_func_get_negative_mu(Anew,S,Si,p_i);

%------------------ IDA (Information Discriminant Analysis)

% disp('Information Discriminant Analysis');

To = Anew;


%To = randn(size(Anew,1),size(Anew,2));
mu_o = 0;

for j = 1:Niter
    %t0 = clock;
    %########### Warning To is not changing throughout iterations
    %Need to update To before each calling of get_megative_mu
    %disp('.');

    %fprintf(1,'.');

    %orthonormalize
    [u s v] = svd(To);
    is = 1./diag(s(1:m,1:m));
    F = 1/(max(max(s)));
    is = diag(is);
    To = is * u' * To;

    mutemp = -dataproc_func_get_negative_mu(To,S,Si,p_i);
    %disp(['initial fcn value: ' num2str(mutemp,'%1.8f')])
    [Topt, mu_opt, exitflag, output, grad] = ...
        fminunc(@(T) dataproc_func_get_negative_mu(T,S,Si,p_i),To,OPT);
    exitflag;
    %disp(['time elapsed: ' num2str(etime(clock,t0))])
    mu_opt = -mu_opt;
    %disp(['fcn value at solution: ' num2str(mu_opt,'%1.8f')])
    %disp(['gradient at solution: ' num2str(norm(grad,inf))])
    if mu_opt > mu_o
        %disp('=====================')
        %disp('better solution found');
        mu_o = mu_opt;
        T_o = Topt;
    end
    To = Topt + randn(size(Topt,1),size(Topt,2));
end

% disp(['optimal function: ' num2str(mu_opt)]);

%orthonormalize
[u s v] = svd(T_o);
is = 1./diag(s(1:m,1:m));
is = diag(is);
Aida = is * u' * T_o;
mu_ida = mu_o;
% 
% figure(2)
% plot(Aida')
% hold on
% pause

% % Features extracted from training data
% % Ftrain {m x sum(Ntr)} = Aida {m x n} * BinData {n x sum(Ntr)}
% Ftrain = Aida * BinData;
% 
% plot(Ftrain(1,find(Memship == 0)), Ftrain(2,find(Memship == 0)), 'r.', Ftrain(1,find(Memship == 1)), Ftrain(2,find(Memship == 1)), 'b.');
% 


Fmat = Aida.';
