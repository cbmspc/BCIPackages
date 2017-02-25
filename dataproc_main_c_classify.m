% function [Decision, Posterior] = icpca_classify (FPmatC,
% DimTrainingDataC, TestData, ClassifierType, Prior, CParm)
%
% Classifier based on Classwise Principle Component Analysis and
% Information Discriminant Analysis
%
% FPmatC = Feature extraction matrix (from high dimensional data to Fdim)
% DimTrainingDataC (cell) = high dimensional training data (dim x obs)
% TestData (vector or matrix) = high dimensional test data (dim x 1)
% ClassifierType = 'knn' for K-nearest-neighbor
%                  'quadratic' to assume normal pdf 
%                  'parzen' to use parzen windows
%
% CParm = K (if using knn), width (if using parzen), not used for quadratic
%
% Prior (vector, 'empirical', or []) = prior probabilities.
%    To set priors manually, Prior = [P(class0) P(class1) P(class2) ...];
%    To set priors to empirical frequencies, Prior = 'empirical';
%    To set priors to equal among all classes, Prior = []; or Prior = '';
% 
% Class label is integer. It starts from 0 and goes positive.
%
% Example:
% Nchan = size(RDS.class0,1);    % number of sensor channels
% Nsamp = size(RDS.class0,2);    % number of digital samples
% Nobs = size(RDS.class0,3);     % number of observatons
% TestData = rand(Nchan,Nsamp);  % use a random matrix as test
% [FPmatC, DimDataC]  = dataproc_main_c_training({RDS.class0, RDS.class1}, 2);
% [Decision, Posterior] = dataproc_main_c_classify(FPmatC, DimDataC, TestData, 'empirical');



function [Decision, Posterior] = dataproc_main_c_classify (FPmatC, DimTrainingDataC, TestData, ClassifierType, Prior, CParm)

%% Gather the dimensions
if iscell(TestData) || length(size(TestData)) < 2
    error('TestData must be a 2D or 3D matrix');
end

if ~iscell(FPmatC) || length(FPmatC) < 2
    error('FPmatC must be a 1D cell of at least length 2.');
end

Nclass = length(DimTrainingDataC);
if Nclass ~= length(FPmatC)
    error('Number of classes indicated by FPmatC and DimTrainingDataC disagree.');
end

for c = 1:Nclass
    NtrialA(c) = size(DimTrainingDataC{c},2);
end


%% Creating or normalizing prior probabilities

if length(who('Prior')) == 0 || length(Prior) == 0 || strcmp(Prior,'equal')  % not specified or left blank
    Prior = ones(1,Nclass);  % Equal prior
elseif strcmp(Prior,'empirical')
    Prior = NtrialA; % Empirical frequencies as prior
end          % Otherwise, user-specified prior
Prior = Prior ./ sum(Prior);  % Normalize

if length(Prior) ~= Nclass
    error('Length of Prior must be equal to Nclass.');
end


%% CParm

if length(who('CParm')) == 0
    CParm = [];
end

%% Reshaping
if length(size(TestData)) == 3
    for t = 1:size(TestData,3)
        DimTestData(:,t) = reshape(TestData(:,:,t), 1, []);
    end
else  % TestData is already reshaped
    DimTestData = TestData;
end
Nobs = size(DimTestData,2);

%% Classifying (classwise)

for s = 1:Nclass
    NormalizationFactor = 0;
    x = (FPmatC{s} * DimTestData ) .' ;

    if strcmp(ClassifierType,'knn')
        Xtrain = [];
        Xlabel = [];
    end

    for c = 1:Nclass
        X = (FPmatC{s} * DimTrainingDataC{c}) .' ;

        if strcmp(ClassifierType,'knn')
            Xlabel = cat(1,Xlabel,(c-1)*ones(size(X,1),1));
            Xtrain = cat(1,Xtrain,X);
        else
            
            switch ClassifierType
                case 'quadratic'
                    f(s,c,:) = dataproc_func_normalpdf(x, X);
                case 'parzen'
                    f(s,c,:) = dataproc_func_parzenwinpdf(x, X, CParm);
                otherwise
                    error(['Unrecognized classifier: ' ClassifierType]);
            end
            
            NormalizationFactor = NormalizationFactor + f(s,c,:) * Prior(c);
        end
    end

    if ~strcmp(ClassifierType,'knn')
        emp1 = find(f(s,1,:)==0);
        emp2 = find(f(s,2,:)==0);
        emp = intersect(emp1,emp2).';
        if length(emp)
            warning(['pdf is zero for both classes at trial #' num2str(emp)]);
        end
    end
    %reshape(NormalizationFactor,1,[])
    
    if strcmp(ClassifierType,'knn')
        [Dec, Pos(:,:,s)] = dataproc_func_knnclassify (x, Xtrain, Xlabel, CParm);
    else
        for c = 1:Nclass
            Posterior(s,c,:) = f(s,c,:) .* Prior(c) ./ NormalizationFactor;
        end
    end
end
clear('NormalizationFactor','Dec','Err');

if strcmp(ClassifierType,'knn')
    Pos = max(Pos,[],3);
    Pos = Pos ./ (sum(Pos,2) * ones(1,Nclass));
    [Y,I] = max(Pos,[],2);
    Dec = I-1;
    Posterior = Pos;
    Decision = Dec;
else
    [Y,I] = max(max(Posterior,[],1));
    Decision = reshape(I-1,[],1);
    clear('Y','I');
    Posterior = reshape(max(Posterior,[],1),Nclass,[])';
    Posterior = Posterior ./ (sum(Posterior,2) * ones(1,Nclass));

end