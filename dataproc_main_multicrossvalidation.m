% function [pcorrect, pconfuse, indivfold, Parameters] =
% dataproc_main_multicrossvalidation (TrainData, TrainLabels, K, M, DRfunC,
% DRparmC, FEfunC, FEparmC, CfunC, CnameC, PriorC, CparmC, SWoutput,
% Trialstotest, SWbalanced)
%
% K-Fold Cross Validation for multiple parameter search
%
% INPUT ARGUMENTS
%
% TrainData: (obs x dim) (Each row is an observation)
%            (Please reshape first if necessary)
%
% TrainLabels: An array of labels, length must == size(TrainData,1)
%
% K: Number of folds. 
%    If K >= 2, uses stratified K number of folds
%    If K <= -2, disables stratification and uses abs(K) as number of folds
%    If abs(K) == 1, calculates in-training classification rate. This will
%     cause overfitting.
%    If K == 0, uses leave-one-out. This takes the longest time to compute.
%
% M: Number of runs. The pie is re-randomized each run.
%
% THE FOLLOWING ARGUMENTS ARE PAIRED CELLS:
% DRfunC and DRparmC must pair (same length)
% FEfunC and FEparmC must pair (same length)
% CfunC, CnameC, and CparmC must pair (same length)
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
% Cname, PriorC, CParm). The output argument must be (Decision) where
% Decision matches one of the elements in TrainLabels.
%
% Cname: The name of classifier you wish to use, passed in as the 4th
% argument into Cfun.
%
% PriorC: Prior probabilities, used in most Bayesian classifiers. This is
% not used in knn, so enter [] for knn. This is also the 5th argument
% passed into the classifier function.
%
% Cparm: Classifier parameter, used by knn to specify number of neighbors
% and by parzen to specify the multiplier of median distance. This is not
% used for classifiers using Gaussian probability density function, so
% enter [] as placeholder. This is also the 6th argument passed into the
% classifier function.
%
% SWoutput: (OPTIONAL, Default = not silent) 
% 0 to silence the progress output (each dot is one fold, each line is one
% run)
%
% Trialstotest: (OPTIONAL, Default = all trials)
% Specify the trial numbers to test. Default (if empty) is all trials.
% Trials are numbered from 1 to size(TrainData,1) in its original order and
% persists through random permutation at the beginning of each run.
% Pcorrect and Pconfuse will only be based on testing these trials.
%
% SWbalanced: (OPTIONAL, Default = 0)
% If set to 1, balance the number of training trials so that all classes
% will have the same number of training trials at all times. Only the first
% N training trials are used, where N is least number of training trials in
% any one class. Because trial sequence is randomized at the beginning of
% each run (except for leave-one-out), different trials will be used for
% training. If set to 2, also balance the number of testing trials.
% 
%
% OUTPUT ARGUMENTS
%
% pcorrect: An array of length M telling the probability of correct
% (weighted by priors). pcorrect = diag(pconfuse) .* Prior
%
% pconfuse: The confusion matrix normalized to probabilities
%
% indivfold: Details of each individual fold. For mth run, kth fold,
% indivfold{dr,fe,cf,pr}{m,k} is a structure with these fields: 'confusion'
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
% 
% Parameters: A structure that stores the parameters used for cross
% validation
%
%
% USAGE EXAMPLE
%
% [pcorrect,pconfuse] = dataproc_main_multicrossvalidation ( ...
% TrainData, TrainLabels, 10, 10, ...
% {'pca','cpca'}, {0.99,0.99}, ...
% {'aida','aida','aida'}, {1,2,3}, ...
% {[],[],[],[]}, {'linear','quadratic','parzen','knn'}, ...
% 'empirical', {0,0,0,3}, 1, [], []);
%
% Here, we do K=10 fold, M=10 run CV. Two dimension reduction
% algorithms (pca first, then cpca (both with parm 0.99)) are performed. 
% Then three feature extraction are performed at different parameters
% (1,2,3) (they are all aida, but can be different if needed). Then four
% classifiers are performed using four corresponding parameters. The
% resulting output variables would be a cell in the format {DR,FE,CF}. 
%
%
% AUTHOR
%
% Po T Wang
% 2008-09-04 (YYYY-MM-DD)
% 2008-10-10 (revision)
% 2009-02-18 adhoc parallel computing feature
% 2010-06-24 fully compatible parallel computing using parfor
% 2011-05-05 parfor within a run
% 2011-07-25 indivconfusion (abandoned)
% 2012-07-31 Add TrialsToTest to non-LOO CV. Removed the restriction to
%            require indivfold for TrialsToTest
% 2012-08-03 Forces TrainLabels to be vertical
% 2012-08-08 Removed the restriction to require K to be <= number of trials
% 2012-08-09 Added SWbalanced
% 2012-08-28 Renamed the input argument variable Opts to SWoutput
%            Fixed a warning message about stratify and balanced
% 2017-04-19 Added SWbalanced for test trials
% 2018-03-08 Updated example usage
%
%
% ptwang@uci.edu
% Biomedical Engineering, University of California Irvine
%


function [pcorrect, pconfuse, indivfold, Parameters] = dataproc_main_multicrossvalidation (TrainData, TrainLabels, K, M, DRfunC, DRparmC, FEfunC, FEparmC, CfunC, CnameC, PriorC, CparmC, SWoutput, Trialstotest, SWbalanced)

%% Getting the parameters

TrainLabels = TrainLabels(:);
Nobs = length(TrainLabels);
if Nobs ~= size(TrainData,1)
    error('Number of observations from TrainData and TrainLabels disagree.');
end

classes = unique(TrainLabels);

Nclass = length(classes);

if Nclass <= 1
    error('At least two classes are required.');
end

NtrialA = zeros(1,Nclass);
for c = 1:Nclass
    NtrialA(c) = length(find(TrainLabels == classes(c)));
end

% By default, all trials are tested once per run.
TestCandidates = 1:sum(NtrialA);


K = round(K);

SWdisablestratify = 0;
if K <= 1 % Disable stratification
    SWdisablestratify = 1;
end

if exist('SWbalanced','var') && ~isempty(SWbalanced) && SWbalanced
else
    SWbalanced = 0;
end

SWleavenothingout = 0;
SWleaveoneout = 0;
if K == 0 % Use leave-one-out
    SWdisablestratify = 1;
    SWleaveoneout = 1;
    %Moved down %K = sum(NtrialA);
elseif abs(K) == 1 % Use leave-nothing-out
    SWleavenothingout = 1;
    SWdisablestratify = 1;
end

M = round(M);
if M < 1
    M = 1;
    warning('DMMCV:MRUN1','Autofixed: M (number of cross validation runs) cannot be less than 1. New value = %i', M);
end

if M > 1
    if SWdisablestratify && K == sum(NtrialA)
        M = 1;
        warning('DMMCV:M1LOOCV','Changed M to 1 due to leave-one-out.');
    elseif SWleavenothingout
        M = 1;
        warning('DMMCV:M1LNOCV','Changed M to 1 due to leave-nothing-out.');
    end
end

DRfunC = internalfunc_forcecell(DRfunC);
DRparmC = internalfunc_forcecell(DRparmC);
FEfunC = internalfunc_forcecell(FEfunC);
FEparmC = internalfunc_forcecell(FEparmC);
CfunC = internalfunc_forcecell(CfunC);
CnameC = internalfunc_forcecell(CnameC);
CparmC = internalfunc_forcecell(CparmC);
PriorC = internalfunc_forcecell(PriorC);

DR = length(DRfunC);
FE = length(FEfunC);
CF = length(CfunC);
PR = length(PriorC);

if DR ~= length(DRparmC)
    error('Lengths of DRfunC and DRparmC differ.');
elseif FE ~= length(FEparmC)
    error('Lengths of FEfunC and FEparmC differ.');
elseif CF ~= length(CnameC) || length(CfunC) ~= length(CparmC)
    error('Lengths of CfunC, CnameC, and/or CparmC differ');
end

for dr = 1:DR
    if strcmp(DRfunC{dr},'cpca')
        DRfunC{dr} = @dataproc_func_cpca;
    elseif strcmp(DRfunC{dr},'pca')
        DRfunC{dr} = @dataproc_func_pca;
    elseif strcmp(DRfunC{dr},'ctpca')
        DRfunC{dr} = @dataproc_func_ctpca;
    elseif strcmp(DRfunC{dr},'fkt')
        DRfunC{dr} = @dataproc_func_fktdr;
    end
end

for fe = 1:FE
    if strcmp(FEfunC{fe},'ida')
        FEfunC{fe} = @dataproc_func_ida;
    elseif strcmp(FEfunC{fe},'aida')
        FEfunC{fe} = @dataproc_func_aida;
    elseif strcmp(FEfunC{fe},'lda')
        FEfunC{fe} = @dataproc_func_lda;
    end
end

for cf = 1:CF
    if isempty(CfunC{cf})
        CfunC{cf} = @dataproc_func_classify;
    end
end

for pr = 1:PR
    if isempty(PriorC{pr})
        PriorC{pr} = 'empirical';
    end
end

SWprogressoutput = 0;
if exist('SWoutput','var') && ~isempty(SWoutput)
    if isnumeric(SWoutput)
        SWprogressoutput = SWoutput;
    else
        SWprogressoutput = 1;
    end
else
    SWoutput = 1;
end

if nargout >= 3
    %|| K == sum(NtrialA)
    SWfindwinnersubspace = 1;
else
    SWfindwinnersubspace = 0;
end

if exist('Trialstotest','var') && ~isempty(Trialstotest)
    %if ~SWfindwinnersubspace
    %    error('To specify which trials to test, indivfold must be requested as output variable.');
    %end
    Trialstotest = unique(round(Trialstotest));
    if ~isempty(find(Trialstotest > sum(NtrialA), 1)) || ~isempty(find(Trialstotest < 1, 1))
        error('Trialstotest contains at least one invalid trial index.');
    end
%     if length(Trialstotest) == length(1:K) && min(Trialstotest == 1:K)
%         warning('DMMCV:ALLTRIALREDUNDANT','Specifying all trials to test is redundant.');
%         SWkoverride = 0;
%     else
%         SWkoverride = 1;
%     end
    TestCandidates = Trialstotest;
    
    if SWprogressoutput
        fprintf('The following %i trials are selected for testing: %s\n', length(TestCandidates), get_contig_groups_string(TestCandidates, ', '));
    end
    
else
    
    if SWprogressoutput
        fprintf('All %i trials are selected for testing. This is the default behavior.\n', length(TestCandidates));
    end
    
    % Default behavior when Trialstotest is not specified
    %SWkoverride = 0;
end

TestCandidates = TestCandidates(:);

NtestcandidateA = zeros(1,Nclass);
for c = 1:Nclass
    NtestcandidateA(c) = nnz(TrainLabels(TestCandidates) == classes(c));
end


if SWleaveoneout
    K = sum(NtestcandidateA);
end

K = abs(K);

if SWdisablestratify && K > sum(NtestcandidateA)
    %K = sum(NtestcandidateA);
    %if K == 1
    %    K = 0;
    %end
    %warning('DMMCV:KFOLDUNSTRAT1','Autofixed: K (number of folds) cannot be more than the total number of testable trials. New value = %i', K);
    warning('DMMCV:KFOLDUNSTRAT1','There are more folds than the number of testable trials.');
elseif ~SWdisablestratify && K > min(NtestcandidateA)
    %K = min(NtestcandidateA);
    %if K == 1
    %    K = 0;
    %end
    %warning('DMMCV:KFOLDSTRAT1','Autofixed: K (number of folds) cannot be more than the minimum number of testable trials among classes. New value = %i', K);
    warning('DMMCV:KFOLDSTRAT1','There are more folds than the number of testable trials per class.');
end

if SWbalanced && SWdisablestratify
    SWbalanced = 0;
    if ~SWleaveoneout && ~SWleavenothingout
        warning('DMMCV:NOBAL', 'Number of trials can only be balanced if stratification is used');
    end
end

if SWprogressoutput
    if SWbalanced == 2
        fprintf('Balanced: For each run, number of training and testing trials per class is equalized.\n');
    elseif SWbalanced == 1
        fprintf('Balanced: For each run, number of training trials per class is equalized.\n');
    end
end

% The variable Trialstotest is re-used for either the list of testable
% trials (leave one out) or the list of folds (non-LOO).
Trialstotest = 1:K;

Parameters.K = K;
Parameters.M = M;
Parameters.DRfunC = DRfunC;
Parameters.DRparmC = DRparmC;
Parameters.FEfunC = FEfunC;
Parameters.FEparmC = FEparmC;
Parameters.CfunC = CfunC;
Parameters.CnameC = CnameC;
Parameters.PriorC = PriorC;
Parameters.CparmC = CparmC;
Parameters.SWoutput = SWoutput;
Parameters.Trialstotest = Trialstotest;
Parameters.TestCandidates = TestCandidates;
Parameters.NtrialA = NtrialA;
Parameters.NtestcandidateA = NtestcandidateA;
Parameters.SWdisablestratify = SWdisablestratify;
Parameters.SWleavenothingout = SWleavenothingout;
Parameters.SWleaveoneout = SWleaveoneout;
Parameters.SWbalanced = SWbalanced;

%% Partitioning data

if SWleavenothingout && SWprogressoutput
    fprintf('Leave-nothing-out (in-training performance) is enabled!\n');
end

NtestA = zeros(1,Nclass);

if SWdisablestratify
    %Ntest = sum(NtrialA)/K;
    Ntest = sum(NtestcandidateA)/K;
    TrainLabelEnums = TrainLabels;
    for c = 1:Nclass
        TrainLabelEnums(TrainLabels==classes(c),1) = (c-1);
    end
    if SWprogressoutput
        fprintf('Stratification is disabled.\n');
    end
    if SWprogressoutput && Ntest == 1
        fprintf('Leave-one-out is enabled.\n');
    elseif SWprogressoutput
        fprintf('Avg # test trials / fold = %g\n',Ntest);
    end
else
    Ntest = [];
    TrainLabelEnums = [];
    idx = cell(1,Nclass);
    for c = 1:Nclass
        % Indices of class 'c'
        idx{c} = find(TrainLabels == classes(c));

        % Number of test (validation) data.
        %NtestA(c) = NtrialA(c)/K;
        NtestA(c) = NtestcandidateA(c)/K;
    end
    if SWprogressoutput
        fprintf('Stratification is enabled.\n');
        fprintf('Avg # test trials / fold = [%s]\n',regexprep(num2str(NtestA), ' +', ', '));
    end
end

confusion = cell(DR,FE,CF,PR);
pconfuse = cell(DR,FE,CF,PR);
pcorrect = cell(DR,FE,CF,PR);
if SWfindwinnersubspace
    indivfold = cell(DR,FE,CF,PR);
end

for dr = 1:DR
    for fe = 1:FE
        for cf = 1:CF
            for pr = 1:PR
                confusion{dr,fe,cf,pr} = zeros(Nclass,Nclass,M);
            end
        end
    end
end

NTrialstotest = length(Trialstotest);

pcmax = zeros(M,1);

for m = 1:M
    if SWprogressoutput
        fprintf('Now doing run #%3i : ',m);
        bksp = 0;
    end
    % Do these every run
    
    if SWdisablestratify
        if Ntest == 1
            % Do leave-one-out sequentially
            ridx = (1:sum(NtrialA)).';
        else
            % Random permute sequence of all labels
            ridx = randperm(sum(NtrialA)).';
        end
    else
        for c = 1:Nclass
            % Random permute the sequence
            r = randperm(NtrialA(c));
            ridx{c} = idx{c}(r);
        end
    end
    %progressmark = 1;
    localconfusion = cell(1,NTrialstotest);
    localindivfold = cell(1,NTrialstotest);
   % %For parfor:
   % %kitick = round(linspace(0,NTrialstotest,21));
   % %kitick = kitick(1:end-1)+1;
    
   % %For non-parallel for:
   % kitick = round(linspace(0,NTrialstotest,21));
   % kitick = kitick(2:end);
   
%    % 20110505
%    if SWprogressoutput
%        kitick = unique(round(linspace(0,NTrialstotest,101)));
%        %kitick = kitick(1:end-1)+1;
%        kiti = 1;
%        spaces = 1+floor(log10(NTrialstotest));
%        InstantProgressFormat = ['%' num2str(spaces) 'i/%' num2str(spaces) 'i'];
%        InstantProgressBackspace = '\b';
%        for di = 1:spaces*2
%            InstantProgressBackspace = [InstantProgressBackspace '\b'];
%        end
%    end
    
    % 75% speed increase using parfor with 4 threads
    %parfor%
    parfor ki = 1:NTrialstotest
        k = Trialstotest(ki); % k is the trial/fold index to be tested
        
        % Do these every fold
        cvtestlabels = [];
        cvtrainlabels = [];
        cvtestdata = [];
        cvtraindata = [];
        %cvtraindatac = [];
        %cvtestdatac = [];
        cvtestdataidx = [];
        cvtraindataidx = [];
        Posterior = [];
        WinnerSubspace = [];
        Tmat = [];
        ridxtest = [];
        ridxtrain = [];
        
        if SWdisablestratify
            [~, IA] = intersect(ridx, TestCandidates);
            tmp = ridx(sort(IA));
            if SWprogressoutput >= 2
                fprintf('RUN %i FOLD %i  TEST SET: %s (%i)\n', m, k, get_contig_groups_string(tmp([round((k-1)*Ntest)+1 : round((k)*Ntest)])', ', '), length(tmp([round((k-1)*Ntest)+1 : round((k)*Ntest)])));
            end
            ridxtest = tmp( round((k-1)*Ntest)+1 : round((k)*Ntest) );
            
            %ridxtest = ridx( round((k-1)*Ntest)+1 : round((k)*Ntest) ); %#ok<PFBNS>
            if SWleavenothingout
                ridxtrain = ridx;
                %ridxtrain = ridxtest;
            else
                ridxtrain = setdiff(ridx,ridxtest);
            end
            if SWprogressoutput >= 2
                fprintf('RUN %i FOLD %i TRAIN SET: %s (%i)\n\n', m, k, get_contig_groups_string(ridxtrain', ', '), length(ridxtrain'));
            end
            
            cvtestlabels = TrainLabelEnums(ridxtest); %#ok<PFBNS>
            cvtrainlabels = TrainLabelEnums(ridxtrain);
            cvtestdata = TrainData(ridxtest,:); %#ok<PFBNS>
            cvtraindata = TrainData(ridxtrain,:);
            cvtestdataidx = ridxtest;
            cvtraindataidx = ridxtrain;
        else
            for c = 1:Nclass
                [~, IA] = intersect(ridx{c}, TestCandidates);
                tmp = ridx{c}(sort(IA));
                if SWprogressoutput >= 2
                    fprintf('RUN %i FOLD %i CLASS %i  TEST SET: %s (%i)\n', m, k, c, get_contig_groups_string(tmp([round((k-1)*NtestA(c))+1 : round((k)*NtestA(c))])', ', '), length(tmp([round((k-1)*NtestA(c))+1 : round((k)*NtestA(c))])));
                end
                
                ridxtest{c} = tmp( round((k-1)*NtestA(c))+1 : round((k)*NtestA(c)) );
                %ridxtest{c} = ridx{c}( round((k-1)*NtestA(c))+1 : round((k)*NtestA(c)) ); %#ok<PFBNS>
                
                ridxtrain{c} = setdiff(ridx{c},ridxtest{c});
            end
            
            %20120809: If SWbalanced, pick least number of trials so
            % training set is balanced among all classes
            if SWbalanced
                tmp = inf;
                for c = 1:Nclass
                    tmp = min(tmp,length(ridxtrain{c}));
                end
                % Pick only the first tmp trials from each class
                for c = 1:Nclass
                    ridxtrain{c} = ridxtrain{c}(1:tmp);
                end
                
                %20170419: If SWbalanced == 2, also balance the test set
                if SWbalanced >= 2
                    tmp = inf;
                    for c = 1:Nclass
                        tmp = min(tmp,length(ridxtest{c}));
                    end
                    % Pick only the first tmp trials from each class
                    for c = 1:Nclass
                        ridxtest{c} = ridxtest{c}(1:tmp);
                    end
                end
            end

            
            for c = 1:Nclass
                if SWprogressoutput >= 2
                    fprintf('RUN %i FOLD %i CLASS %i TRAIN SET: %s (%i)\n\n', m, k, c, get_contig_groups_string(ridxtrain{c}', ', '), length(ridxtrain{c}'));
                end
                
                cvtestlabels = cat(1,cvtestlabels,(c-1)*ones(length(ridxtest{c}),1));
                cvtrainlabels = cat(1,cvtrainlabels,(c-1)*ones(length(ridxtrain{c}),1));
                cvtestdata = cat(1,cvtestdata,TrainData(ridxtest{c},:));
                cvtraindata = cat(1,cvtraindata,TrainData(ridxtrain{c},:));
                cvtestdataidx = cat(1,cvtestdataidx,ridxtest{c});
                cvtraindataidx = cat(1,cvtraindataidx,ridxtrain{c});
            end
        end
        
        for dr = 1:DR
            cvtraindatac = [];
            cvtestdatac = [];
            DRfun = DRfunC{dr}; %#ok<PFBNS>
            DRparm = DRparmC{dr}; %#ok<PFBNS>

            % Dimension reduction
            if ~isempty(DRfun)
                % DRmat can potentially be a cell. Always treat it as a cell
                DRmat = DRfun(cvtraindata,cvtrainlabels,DRparm);
                if ~iscell(DRmat)
                    %tmp = DRmat;
                    %clear('DRmat');
                    DRmat = {DRmat};
                    %clear('tmp');
                end
                if ~isempty(DRmat)
                    S = length(DRmat);
                    cvtraindatac = cell(1,S);
                    cvtestdatac = cell(1,S);
                    for s = 1:S
                        cvtraindatac{s} = cvtraindata * DRmat{s};
                        cvtestdatac{s} = cvtestdata * DRmat{s};
                    end
                end
            else
                DRmat = {};
                cvtraindatac = {cvtraindata};
                cvtestdatac = {cvtestdata};
            end
            
            
            
            
            for fe = 1:FE
                cvftraindatac = cell(1);
                cvftestdatac = cell(1);
                FEfun = FEfunC{fe}; %#ok<PFBNS>
                FEparm = FEparmC{fe}; %#ok<PFBNS>
               
                % Feature extraction
                if ~isempty(FEfun)
                    S = length(DRmat);
                    FEmat = {FEfun(cvtraindatac{1},cvtrainlabels,FEparm)};
                    cvftraindatac{1} = cvtraindatac{1} * FEmat{1};
                    cvftestdatac{1} = cvtestdatac{1} * FEmat{1};
                    for s = 2:S
                        FEmat{s} = FEfun(cvtraindatac{s},cvtrainlabels,FEparm);
                        cvftraindatac{s} = cvtraindatac{s} * FEmat{s};
                        cvftestdatac{s} = cvtestdatac{s} * FEmat{s};
                    end
                else
                    %FEmat = {};
                    cvftraindatac = cvtraindatac;
                    cvftestdatac = cvtestdatac;
                end


                for cf = 1:CF
                    Cfun = CfunC{cf}; %#ok<PFBNS>
                    Cname = CnameC{cf}; %#ok<PFBNS>
                    Cparm = CparmC{cf}; %#ok<PFBNS>
                    for pr = 1:PR
                        Prior = PriorC{pr}; %#ok<PFBNS>
                        
                        % Classification
                        
                        if SWfindwinnersubspace
                            %[Decision,Posterior,WinnerSubspace] = Cfun(cvftestdatac,cvftraindatac,cvtrainlabels,Cname,Prior,Cparm);
                            [Decision,Posterior,WinnerSubspace] = Cfun(cvftestdatac,cvftraindatac,cvtrainlabels,Cname,Prior,Cparm);
                        else
                            %[size(cvftestdatac{1}) size(cvftestdatac{2}) size(cvftraindatac{1}) size(cvftraindatac{2})]
                            Decision = Cfun(cvftestdatac,cvftraindatac,cvtrainlabels,Cname,Prior,Cparm);
                        end
                        
                        if size(Decision,2) > size(Decision,1)
                            Decision = Decision.';
                        end
                        
                        %if SWfindwinnersubspace
                        %    correctall = [];
                        %end
                        
                        for x = 0:Nclass-1
                            for y = 0:Nclass-1
                                %intersection = intersect(find(Decision == x), find(cvtestlabels == y));
                                Tmat.confusion(x+1,y+1) = length(intersect(find(Decision == x), find(cvtestlabels == y)));
                                %20100623% confusion{dr,fe,cf,pr}(x+1,y+1,m) = confusion{dr,fe,cf,pr}(x+1,y+1,m) + Tmat.confusion(x+1,y+1);
                                %if x == y & SWfindwinnersubspace
                                %    correctall = union(correctall,intersection);
                                %end
                            end
                        end
                        
                        localconfusion{ki}{dr,fe,cf,pr} = Tmat.confusion;
                        
                        if SWfindwinnersubspace
                            Tmat.traintable = [cvtraindataidx cvtrainlabels];
                            Tmat.testtable = [cvtestdataidx cvtestlabels Decision WinnerSubspace-1 Posterior];
                            %Tmat.trfmat{1} = eye(size(cvtestdata,2));
                            %if length(DRmat) > 0
                            %    Tmat.trfmat{1} = Tmat.trfmat{1} * DRmat{1};
                            %end
                            %if length(FEmat) > 0
                            %    Tmat.trfmat{1} = Tmat.trfmat{1} * FEmat{1};
                            %end
                            %for s = 2:length(DRmat)
                            %    Tmat.trfmat{s} = eye(size(cvtestdata,2));
                            %    Tmat.trfmat{s} = Tmat.trfmat{s} * DRmat{s};
                            %    if length(FEmat) > 0
                            %        Tmat.trfmat{s} = Tmat.trfmat{s} * FEmat{s};
                            %    end
                            %end
                            %20100623% indivfold{dr,fe,cf,pr}{m,ki} = Tmat;
                            localindivfold{ki}{dr,fe,cf,pr} = Tmat;
                        end
                    end % End of a prior
                end % End of a classifier
            end % End of a feature extraction
        end % End of a dimension reduction

        %while SWprogressoutput && ki/NTrialstotest*10 >= progressmark
        %    fprintf('.');
        %    progressmark = progressmark + 1;
        %end
        
        %if SWprogressoutput && ~isempty(intersect(ki,kitick))
        %    fprintf('.');
        %end
        
%         %20110505
%         if SWprogressoutput
%             if ki >= kitick(kiti)
%                 if bksp
%                     fprintf(InstantProgressBackspace);
%                 end
%                 bksp = 1;
%                 fprintf(InstantProgressFormat, ki, NTrialstotest);
%                 kiti = kiti+1;
%             end
%         end
        
        %For dots:
        %if SWprogressoutput && ~isempty(intersect(ki,kitick))
        %    for di = 1:length(find(ki==kitick))
        %        fprintf('.');
        %    end
        %end
        
    end % End of a fold (ki)
    
    % Collect data after parallel loop
    %20100623% confusion{dr,fe,cf,pr}(x+1,y+1,m) = confusion{dr,fe,cf,pr}(x+1,y+1,m) + Tmat.confusion(x+1,y+1);
    %localconfusion{ki}{dr,fe,cf,pr}{m} = Tmat.confusion;
    for dr = 1:DR
        for fe = 1:FE
            for cf = 1:CF
                for pr = 1:PR
                    for ki = 1:NTrialstotest
                        confusion{dr,fe,cf,pr}(:,:,m) = confusion{dr,fe,cf,pr}(:,:,m) + localconfusion{ki}{dr,fe,cf,pr};
                    end
                end
            end
        end
    end
    
    %20100623% indivfold{dr,fe,cf,pr}{m,ki} = Tmat;
    %localindivfold{ki}{dr,fe,cf,pr}{m} = Tmat;
    if SWfindwinnersubspace
        for dr = 1:DR
            for fe = 1:FE
                for cf = 1:CF
                    for pr = 1:PR
                        for ki = 1:NTrialstotest
                            indivfold{dr,fe,cf,pr}{m,ki} = localindivfold{ki}{dr,fe,cf,pr};
                        end
                    end
                end
            end
        end
    end
    for dr = 1:DR
        for fe = 1:FE
            for cf = 1:CF
                for pr = 1:PR
                    %Prior = PriorC(pr);
                    pconfuse{dr,fe,cf,pr}(:,:,m) = confusion{dr,fe,cf,pr}(:,:,m) ./ (ones(size(confusion{dr,fe,cf,pr}(:,:,m),1),1) * sum(confusion{dr,fe,cf,pr}(:,:,m),1));
                    Nactualclasstrial = sum(confusion{dr,fe,cf,pr}(:,:,m),1);
                    %Ntotal = sum(Nactualclasstrial);
                    %pcorrect{dr,fe,cf,pr}(m) = 0;
                    %for i = 1:Nclass
                    %    %pcorrect{dr,fe,cf,pr}(m) = pcorrect{dr,fe,cf,pr}(m) + confusion{dr,fe,cf,pr}(i,i,m) / Nactualclasstrial(i) * Nactualclasstrial(i)/Ntotal;
                    %end
                    %20100812 %pcorrect{dr,fe,cf,pr}(m) = sum(diag(pconfuse{dr,fe,cf,pr}(:,:,m)).' .* prior);
                    % A more direct pcorrect calculation
                    pcorrect{dr,fe,cf,pr}(m) = sum(diag(confusion{dr,fe,cf,pr}(:,:,m))) ./ sum(sum(confusion{dr,fe,cf,pr}(:,:,m)),2);
                end % End of PR
            end % End of CF
        end % End of FE
    end % End of DR
    if SWprogressoutput % && ~SWkoverride
        ChancePcorrect = NaN;
        for dr = 1:DR
            for fe = 1:FE
                for cf = 1:CF
                    for pr = 1:PR
                        Prior = PriorC{pr};
                        if pcorrect{dr,fe,cf,pr}(m) > pcmax(m)
                            pcmax(m) = pcorrect{dr,fe,cf,pr}(m);
                            if strcmp(Prior,'empirical')
                                prior = Nactualclasstrial ./ sum(Nactualclasstrial);
                            elseif strcmp(Prior,'equal')
                                prior = ones(1,Nclass) ./ Nclass;
                            else
                                prior = Prior ./ sum(Prior);
                            end
                            ChancePcorrect = max(prior);
                        end
                    end
                end
            end
        end
        if bksp
            fprintf([InstantProgressBackspace '\b']); %#ok<UNRCH>
        end
        fprintf(' Best Pcorrect = %-5.3f (chance = %-5.3f)\n',pcmax(m),ChancePcorrect);
    end
end % End of a run (m)

% if SWkoverride
%     pcorrect = {[]};
%     pconfuse = {[]};
% end

% 2011-07-25: There is no "variance" in leave-one-out, because the
% procedure is deterministic
%
% if SWfindwinnersubspace
%     indivconfusion = cell(DR,FE,CF,PR);
%     indivpconfuse = cell(DR,FE,CF,PR);
%     indivpcorrect = cell(DR,FE,CF,PR);
%     for dr = 1:DR
%         for fe = 1:FE
%             for cf = 1:CF
%                 for pr = 1:PR
%                     kmax = numel(indivfold{dr,fe,cf,pr});
%                     indivconfusion{dr,fe,cf,pr} = zeros(Nclass,Nclass,kmax);
%                     indivpconfuse{dr,fe,cf,pr} = zeros(Nclass,Nclass,kmax);
%                     indivpcorrect{dr,fe,cf,pr} = zeros(1,kmax);
%                     for k = 1:kmax
%                         indivconfusion{dr,fe,cf,pr}(:,:,k) = indivfold{dr,fe,cf,pr}{k}.confusion;
%                     end
%                     divisor = ones(Nclass,1)*sum(mean(indivconfusion{dr,fe,cf,pr},3),1);
%                     for k = 1:kmax
%                         indivpconfuse{dr,fe,cf,pr}(:,:,k) = indivconfusion{dr,fe,cf,pr}(:,:,k) ./ divisor;
%                         indivpcorrect{dr,fe,cf,pr}(k) = sum(diag(indivconfusion{dr,fe,cf,pr}(:,:,k))) ./ sum(sum(indivconfusion{dr,fe,cf,pr}(:,:,k),1),2);
%                     end
%                 end
%             end
%         end
%     end
%     % standard deviation does not make sense?
%     pconfuse = indivpconfuse;
%     pcorrect = indivpcorrect;
% end



function V = internalfunc_forcecell (V)
if ~iscell(V)
    V = {V};
%     tmp = V;
%     clear('V');
%     V{1} = tmp;
end

