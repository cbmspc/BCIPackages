% function [pcorrect, pconfuse, indivfold] =
% dataproc_main_crossvalidation (TrainData, TrainLabels, K, M, DRfun,
% DRparm, FEfun, FEparm, Cfun, Cname, Prior, Cparm)
%
% Stratified K-Fold Cross Validation
%
% INPUT ARGUMENTS
%
% TrainData: (obs x dim) (Each row is an observation)
%            (Please reshape first if necessary)
%
% TrainLabels: An array, length must == size(TrainData,1)
%
% K: Number of folds. K=0 to use leave-one-out. K<-1 to disable
% stratification.
%
% M: Number of runs. The pie is re-randomized each run.
%
% DRfun: Dimension reduction function. If left blank [] or '', DR will not
% be performed. Input arguments must be (TrainData, TrainLabels, DRparm).
% Output argument must be a transformation matrix T such that
% TrainData * T = ReducedDimensionTrainData. Output argument can also be a
% cell, in which each element will be treated as a subspace transformation.
% Feature extraction will be performed as many times as the number of
% subspaces. NB: If time-bin or whitening is required, do them before using
% this cross validation function.
%
% DRparm: The parameter for DRfun, passed as the third argument
%
% FEfun: Feature extraction function. If left blank [] or '', FE will not
% be performed. Input arguments must be (TrainData, TrainLabels, FEparm).
% TrainData can be the original high-dimensional TrainData or the reduced
% dimension data. Output argument must be a transformation matrix T such
% that TrainData * T = FeatureTrainData
%
% Cfun: Classifier function. If left blank [] or '', the default classifier
% will be used. Input arguments must be (TestData, TrainData, TrainLabels,
% Cname, Prior, CParm). The output argument must be (Decision) where
% Decision matches one of the elements in TrainLabels.
%
% Cname: The name of classifier you wish to use, passed in as the 4th
% argument into Cfun.
%
% Prior: Prior probabilities, used in most Bayesian classifiers. This is
% not used in knn, so enter [] for knn. This is also the 5th argument
% passed into the classifier function.
%
% Cparm: Classifier parameter, used by knn to specify number of neighbors
% and by parzen to specify the multiplier of median distance. This is not
% used for classifiers using Gaussian probability density function, so
% enter [] as placeholder. This is also the 6th argument passed into the
% classifier function.
%
% Opts: 0 to silence the progress output (each dot is one fold, each line
% is one run)
%
% OUTPUT ARGUMENTS
%
% pcorrect: An array of length M telling the probability of correct.
%
% pconfuse: The confusion matrix normalized to probabilities
%
% indivfold: Details of each individual fold. For mth run, kth fold,
% indivfold{dr,fe,cf}{m,k} is a structure with these fields: 'confusion'
% is the confusion matrix of this fold. 'traintable' is a list of (1)
% indices and (2) their label indices (starts with 0) used for training
% this fold. 'testtable' is a table of four plus Nclass column vectors:
% (1) indices used for testing (correspond to the original sequence in
% TrainData and TrainLabels, (2) truth label index (starts with 0),
% (3) decision label index (starts with 0) by classifier, (4) index
% (starts with 0) of the subspace that was used to reach the decision.
% (5...) Posterior probabilities
%
%
% AUTHOR
%
% Po T Wang
% 2008-08-29 (YYYY-MM-DD)
% 2008-10-10 (revision)
% ptwang@uci.edu
% Biomedical Engineering, University of California Irvine
%


function [pcorrect, pconfuse, indivfold] = crossvalidation (TrainData, TrainLabels, K, M, DRfun, DRparm, FEfun, FEparm, Cfun, Cname, Prior, Cparm, Opts, Koverride)

if length(who('Opts')) == 0
    Opts = 1;
end

if isempty(who('Koverride'))
    Koverride = [];
end

SWsilent = 0;
if Opts == -1
    SWsilent = 1;
    Opts = 0;
end

if nargout >= 3
    [pcorrect,pconfuse,indivfold] = dataproc_main_multicrossvalidation(TrainData,TrainLabels,K,M,DRfun,DRparm,FEfun,FEparm,Cfun,Cname,Prior,Cparm,Opts,Koverride);
    indivfold = indivfold{1,1,1};
else
    [pcorrect,pconfuse] = dataproc_main_multicrossvalidation(TrainData,TrainLabels,K,M,DRfun,DRparm,FEfun,FEparm,Cfun,Cname,Prior,Cparm,Opts);
end

pcorrect = pcorrect{1,1,1};
pconfuse = pconfuse{1,1,1};

classes = unique(TrainLabels);
Nclass = length(classes);
for c = 1:Nclass
    Ntrial(c) = length(find(TrainLabels==classes(c)));
end

if SWsilent == 0
    disp(' ');
    disp(' ');


    Kfoldcomment = '';
    if K == 0
        Kfoldcomment = '(leave-one-out)';
    elseif abs(K) == 1
        Kfoldcomment = '(in-training-rate, not CV!)';
    elseif K < -1
        Kfoldcomment = '(unstratified CV)';
    elseif K > 1
        Kfoldcomment = '(stratified CV)';
    end

    disp('Cross Validation Report');
    disp(['   # runs * # folds: ' num2str(M) ' * ' num2str(K) ' ' Kfoldcomment]);
    disp(['     Data dimension: ' num2str(size(TrainData,2))]);
    disp(['# samples per class: ' num2str(Ntrial)]);
    disp(['Dimension reduction: ' DRfun '(' num2str(DRparm) ')']);
    disp([' Feature extraction: ' FEfun '(' num2str(FEparm) ')']);
    if length(Cfun) == 0
        disp(['Classifier function: ' 'Bayesian']);
    else
        disp(['Classifier function: ' Cfun]);
    end
    disp(['   Density function: ' Cname '(' num2str(Cparm) ')']);
    disp(['Prior probabilities: ' '[' num2str(Prior) ']']);
    disp(' ');
    disp(['Classification rate: ' num2str(mean(pcorrect*100),3) '% ± ' num2str(std(pcorrect*100),3) '%']);
    disp(' ');
    disp('Confusion matrix');
    
    %disp([DRfun '(' num2str(DRparm) ')' '->' FEfun '(' num2str(FEparm) ')' '->' Cfun '/' Cname '(' num2str(Cparm) ',' '[' num2str(Prior) ']' ')'])

    me = mean(pconfuse,3) * 100;
    sd = std(pconfuse,[],3) * 100;
    
    for i = 1:size(me,1)
        for j = 1:size(me,2)
            fprintf('%.1f%% ± %.1f%% \t', me(i,j), sd(i,j));
        end
        fprintf('\n');
    end
    
    
    disp(' ');
    disp(' ');
end

% %% Getting the parameters
% Nobs = length(TrainLabels);
% if Nobs ~= size(TrainData,1)
%     error('Number of observations from TrainData and TrainLabels disagree.');
% end
%
% classes = unique(TrainLabels);
%
% Nclass = length(classes);
%
% if Nclass <= 1
%     error('At least two classes are required.');
% end
%
% for c = 1:Nclass
%     NtrialA(c) = length(find(TrainLabels == classes(c)));
% end
%
%
% K = round(K);
% if K < 2
%     error('K (number of folds) cannot be less than 2');
% end
%
% M = round(M);
% if M < 1
%     error('M (number of cross validation runs) cannot be less than 1');
% end
%
% if strcmp(DRfun,'cpca')
%     DRfun = @dataproc_func_cpca;
% elseif strcmp(DRfun,'pca')
%     DRfun = @dataproc_func_pca;
% elseif strcmp(DRfun,'ctpca')
%     DRfun = @dataproc_func_ctpca;
% end
%
% if strcmp(FEfun,'ida')
%     FEfun = @dataproc_func_ida;
% elseif strcmp(FEfun,'aida')
%     FEfun = @dataproc_func_aida;
% elseif strcmp(FEfun,'lda')
%     FEfun = @dataproc_func_lda;
% end
%
%
% if length(Cfun) == 0
%     Cfun = @dataproc_func_classify;
% end
%
% % if strcmp(Cfun,'knn') || strcmp(Cfun,'parzen') || ...
% %         strcmp(Cfun,'quadratic') || strcmp(Cfun,'linear') || ...
% %         strcmp(Cfun,'diagquadratic') || strcmp(Cfun,'diaglinear')
% %     Cfun = 'dataproc_func_classify';
% % end
%
%
% %% Partitioning data
% for c = 1:Nclass
%     % Indices of class 'c'
%     idx{c} = find(TrainLabels == classes(c));
%
%     % Number of test (validation) data.
%     NtestA(c) = ceil(NtrialA(c)/K);
% end
%
% confusion = zeros(Nclass,Nclass,M);
%
% for m = 1:M
%     % Do these every run
%     for c = 1:Nclass
%         % Random permute the sequence
%         r = randperm(NtrialA(c));
%         ridx{c} = idx{c}(r);
%     end
%
%     for k = 1:K
%         % Do these every fold
%         cvtestlabels = [];
%         cvtrainlabels = [];
%         cvtestdata = [];
%         cvtraindata = [];
%         cvtraindatac = [];
%         cvtestdatac = [];
%
%         for c = 1:Nclass
%             ridxtest{c} = ridx{c}( (k-1)*NtestA(c)+1 : min([k*NtestA(c),NtrialA(c)]) );
%             ridxtrain{c} = setdiff(ridx{c},ridxtest{c});
%             cvtestlabels = cat(1,cvtestlabels,(c-1)*ones(length(ridxtest{c}),1));
%             cvtrainlabels = cat(1,cvtrainlabels,(c-1)*ones(length(ridxtrain{c}),1));
%             cvtestdata = cat(1,cvtestdata,TrainData(ridxtest{c},:));
%             cvtraindata = cat(1,cvtraindata,TrainData(ridxtrain{c},:));
%         end
%
%         % Dimension reduction
%         if length(DRfun)
%             % DRmat can potentially be a cell. Always treat it as a cell
%             DRmat = DRfun(cvtraindata,cvtrainlabels,DRparm);
%             if ~iscell(DRmat)
%                 tmp = DRmat;
%                 clear('DRmat');
%                 DRmat{1} = tmp;
%                 clear('tmp');
%             end
%             for s = 1:length(DRmat)
%                 cvtraindatac{s} = cvtraindata * DRmat{s};
%                 cvtestdatac{s} = cvtestdata * DRmat{s};
%             end
%         else
%             DRmat = {};
%             cvtraindatac{1} = cvtraindata;
%             cvtestdatac{1} = cvtestdata;
%         end
%
%         % Feature extraction
%         if length(FEfun)
%             FEmat{1} = FEfun(cvtraindatac{1},cvtrainlabels,FEparm);
%             cvtraindatac{1} = cvtraindatac{1} * FEmat{1};
%             cvtestdatac{1} = cvtestdatac{1} * FEmat{1};
%             for s = 2:length(DRmat)
%                 FEmat{s} = FEfun(cvtraindatac{s},cvtrainlabels,FEparm);
%                 cvtraindatac{s} = cvtraindatac{s} * FEmat{s};
%                 cvtestdatac{s} = cvtestdatac{s} * FEmat{s};
%             end
%         else
%             FEmat = {};
%         end
%
%         % Classification
%         Decision = Cfun(cvtestdatac,cvtraindatac,cvtrainlabels,Cname,Prior,Cparm);
%
%         if size(Decision,2) > size(Decision,1)
%             Decision = Decision.';
%         end
%
%         for x = 0:Nclass-1
%             for y = 0:Nclass-1
%                 confusion(x+1,y+1,m) = confusion(x+1,y+1,m) + length(intersect(find(Decision == x), find(cvtestlabels == y)));
%             end
%         end
%     end % End of a fold
%     pconfuse(:,:,m) = confusion(:,:,m) ./ (ones(size(confusion(:,:,m),1),1) * sum(confusion(:,:,m),1));
%     Nactualclasstrial = sum(confusion(:,:,m),1);
%     Ntotal = sum(Nactualclasstrial);
%     pcorrect(m) = 0;
%     for i = 1:Nclass
%         pcorrect(m) = pcorrect(m) + confusion(i,i,m) / Nactualclasstrial(i) * Nactualclasstrial(i)/Ntotal;
%     end
% end % End of a run
%
