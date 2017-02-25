% Economical Locally Linear Embedding, a dimension reduction algorithm
% originally by Sam T Roweis and Lawrence K Saul. For regression.
% This version is optimized for computation and real-time decoding.
%
% To TRAIN the function, blah
%
% To REDUCE a new data ("test data"), blah
%
%
%
% All data are stored as N by D: 
%   N = number of observations
%   D = number of dimensions
%
% Input arguments:
%   Xtrain = training data
%   Xtest  = test (new) data
%   K = K nearest G
%   L = Reduced number of dimensions (for output)
%   Rtrain = The "ground truth" signal. It must have N rows like Xtrain.
%      If this is specified, neighbor search will use it.
%      Otherwise, neighbors will be searched from the training space.
%
% Output arguments:
%   Y = data projected in reduced dimension
%     NOTE: If sufficiently accurate estimate cannot be obtained, Y will be
%           empty. You should increase K and try again.
%   W = weight matrix
%  NS = Nearest G search structure
%

function [Y, TrainObjStruct] = dataproc_func_elle (Xtrain, K, L, Xtest, TrainObjStruct, Rtrain, SWapprox)

if ~exist('SWapprox','var') || isempty(SWapprox)
    SWapprox = 1;
end

if ~exist('TrainObjStruct','var') || isempty(TrainObjStruct) || ~isstruct(TrainObjStruct)
    TrainObjStruct = struct;
end

N = size(Xtrain,1);
D = size(Xtrain,2);

if ~exist('Xtest','var') || isempty(Xtest)
    Xtest = [];
end
Ntest = size(Xtest,1);

if ~exist('Rtrain','var') || isempty(Rtrain)
    Rtrain = [];
end
R = size(Rtrain,1);


%% Steps 1 & 2 TRAIN -- Find K nearest neighbors and the distances
%                    -- Compute best reconstruction weights

% This is required because it is used for decoding new data
if isfield(TrainObjStruct,'NS')
    NS = TrainObjStruct.NS;
elseif D < 10
    NS = KDTreeSearcher(Xtrain);
else
    NS = ExhaustiveSearcher(Xtrain);
end

% If Rtrain is specified, theta space is used to locate neighbors during
% training
if R
    if isfield(TrainObjStruct,'RNS')
        RNS = TrainObjStruct.RNS;
    elseif size(Rtrain,2) < 10
        RNS = KDTreeSearcher(Rtrain);
    else
        RNS = ExhaustiveSearcher(Rtrain);
    end
end

if isfield(TrainObjStruct,'Gtrain') && size(TrainObjStruct.Gtrain,1) == N && size(TrainObjStruct.Gtrain,2) == K
    Gtrain = TrainObjStruct.Gtrain;
else
    if R
        Gtrain = knnsearch(RNS,Rtrain,'k',K+1);
        Gtrain = Gtrain(:,2:end);
    else
        Gtrain = knnsearch(NS,Xtrain,'k',K+1);
        Gtrain = Gtrain(:,2:end);
    end
end

tol = 0;
if K > D
    % Regularization is necessary if K > D
    tol = 1e-3;
end


if isfield(TrainObjStruct,'Wtrain') && size(TrainObjStruct.Wtrain,1) == N && size(TrainObjStruct.Wtrain,2) == K
    Wtrain = TrainObjStruct.Wtrain;
else
    Wtrain = zeros(N,K);
    for n=1:N
        z = Xtrain(Gtrain(n,:),:) - ones(K,1)*Xtrain(n,:);
        c = z * z';
        if tol > 0
            c = c + eye(K,K)*tol*trace(c);
        end
        Wtrain(n,:) = c \ ones(K,1);
        Wtrain(n,:) = Wtrain(n,:) / sum(Wtrain(n,:));
    end
end

%% Steps 1 & 2 TEST  -- Find K nearest neighbors and the distances
%                    -- Compute best reconstruction weights

if ~SWapprox && Ntest
    Gtest = knnsearch(NS,Xtest,'k',K);
    Wtest = zeros(Ntest,K);
    for n=1:Ntest
        z = Xtrain(Gtest(n,:),:) - ones(K,1)*Xtest(n,:);
        c = z * z';
        if tol > 0
            c = c + eye(K,K)*tol*trace(c);
        end
        Wtest(n,:) = c \ ones(K,1);
        Wtest(n,:) = Wtest(n,:) / sum(Wtest(n,:));
    end
else
    Gtest = [];
    Wtest = [];
end

%% Step 3 -- Compute embedding from eigenvectors of cost matrix
if Ntest && SWapprox
    if isfield(TrainObjStruct,'Ytrain') && size(TrainObjStruct.Ytrain,1) == N && size(TrainObjStruct.Ytrain,2) == L
        Ytrain = TrainObjStruct.Ytrain;
    else
        % Ytrain was not calculated.
        M = Build_M (TrainObjStruct, N, K, 0, Wtrain, Gtrain, [], []);
        Ytrain = Build_Y (M, L, N, 0);
    end
    Y = Build_Y_approx (Ytrain, NS, Xtrain, Xtest, Ntest, L, D);
else
    M = Build_M (TrainObjStruct, N, K, Ntest, Wtrain, Gtrain, Wtest, Gtest);
    Y = Build_Y (M, L, N, Ntest);
end

%% Save data
TrainObjStruct.Wtrain = Wtrain;
TrainObjStruct.NS = NS;
TrainObjStruct.Gtrain = Gtrain;
if ~Ntest
    TrainObjStruct.Mtrain = M;
    TrainObjStruct.Ytrain = Y;
end


function M = Build_M (TrainObjStruct, N, K, Ntest, Wtrain, Gtrain, Wtest, Gtest)
%% Builds M, to be decomposed in the next step
if isfield(TrainObjStruct,'Mtrain') && size(TrainObjStruct.Mtrain,1) == size(TrainObjStruct.Mtrain,2) && size(TrainObjStruct.Mtrain,1) == N
    % We already have M for training data. We only need to update it with
    % test data. Prepare to build M from N+1
    M = sparse(1:N+Ntest,1:N+Ntest,ones(1,N+Ntest),N+Ntest,N+Ntest,4*K*(N+Ntest));
    M(1:N,1:N) = TrainObjStruct.Mtrain;
    nstart = N+1;
else
    % Prepare to build M from scratch
    M = sparse(1:N+Ntest,1:N+Ntest,ones(1,N+Ntest),N+Ntest,N+Ntest,4*K*(N+Ntest));
    nstart = 1;
end

for n=nstart:N+Ntest
    if n>N
        w = Wtest(n-N,:);
        m = Gtest(n-N,:);
    else
        w = Wtrain(n,:);
        m = Gtrain(n,:);
    end
    M(n,m) = M(n,m) - w; %#ok<*SPRIX>
    M(m,n) = M(m,n) - w';
    M(m,m) = M(m,m) + w'*w;
end


function Y = Build_Y (M, L, N, Ntest)
%% Builds embedding from scratch
opts.disp = 0;
opts.isreal = 1;
opts.issym = 1;
[Y v] = eigs(M, L+1, 0, opts);
v = diag(v);
[v i] = sort(v, 'descend');
Y = Y(:,i);
if v(end) > eps
    %warning('ELLE:NOZERO','Data is too noisy. Resulting embedding may be too inaccurate.');
end
Y = Y(:,1:end-1);

for i = 1:size(Y,2)
    if sign(Y(find(Y(:,i),1),i)) < 1
        Y(:,i) = -Y(:,i); % force the first non-zero element to be positive
    end
end
if Ntest
    Y = Y(N+1:end,:);
end
%Y = Y * sqrt(N+Ntest);



function Ytest = Build_Y_approx (Ytrain, NS, Xtrain, Xtest, Ntest, L, D)
Method = 2;
if Method == 1
    %% Nearest neighbor method
    Gtest = knnsearch(NS,Xtest,'k', 1);
    Ytest = Ytrain(Gtest,:);
    
elseif Method == 2
    %% Triangle approximation method
    Gtest = knnsearch(NS,Xtest,'k', L+1);
    Ytest = zeros(Ntest, L);
    for n = 1:Ntest
        % Create basis vectors
        v = zeros(D,L+1);
        for k = 1:L+1
            v(:,k) = Xtrain(Gtest(n,k),:).';
        end
        vv = zeros(D,L);
        % Re-reference to the first neighbor point
        for k = 2:L+1
            vv(:,k-1) = v(:,k) - v(:,1);
        end
        alpha = vv \ (Xtest(n,:).' - v(:,1));
        u = zeros(L,L+1);
        for k = 1:L+1
            u(:,k) = Ytrain(Gtest(n,k),:).';
        end
        uu = zeros(L,L);
        for k = 2:L+1
            % Re-reference to the first neighbor point
            uu(:,k-1) = u(:,k) - u(:,1);
        end
        Ytest(n,:) = uu*alpha + u(:,1);
    end
end
