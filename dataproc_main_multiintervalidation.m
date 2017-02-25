% function [pcorrect, pconfuse, indivfold] =
% dataproc_main_multiintervalidation (Data, Labels, TestIdx, DRfunC,
% DRparmC, FEfunC, FEparmC, CfunC, CnameC, Prior, CparmC, Opts)
%
% Inter Validation for multiple parameter search
%
% INPUT ARGUMENTS
%
% Data: (obs x dim) (Each row is an observation)
%            (Please reshape first if necessary)
%
% Labels: An array of labels, length must == size(Data,1)
%
% THE FOLLOWING ARGUMENTS ARE PAIRED CELLS:
% DRfunC and DRparmC must pair (same length)
% FEfunC and FEparmC must pair (same length)
% CfunC, CnameC, and CparmC must pair (same length)
%
% DRfun: Dimension reduction function. If left blank [] or '', DR will not
% be performed. Input arguments must be (Data, Labels, DRparm).
% Output argument must be a transformation matrix T such that
% Data * T = ReducedDimensionData. Output argument can also be a
% cell, in which each element will be treated as a subspace transformation.
% Feature extraction will be performed as many times as the number of
% subspaces. NB: If time-bin or whitening is required, do them before using
% this cross validation function.
%
% DRparm: The parameter for DRfun, passed as the third argument
%
% FEfun: Feature extraction function. If left blank [] or '', FE will not
% be performed. Input arguments must be (Data, Labels, FEparm).
% Data can be the original high-dimensional Data or the reduced
% dimension data. Output argument must be a transformation matrix T such
% that Data * T = FeatureData
%
% Cfun: Classifier function. If left blank [] or '', the default classifier
% will be used. Input arguments must be (TestData, Data, Labels,
% Cname, Prior, CParm). The output argument must be (Decision) where
% Decision matches one of the elements in Labels.
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
%
% OUTPUT ARGUMENTS
%
% pcorrect: An array of length M telling the probability of correct
% (weighted by priors). pcorrect = diag(pconfuse) .* Prior
%
% pconfuse: The confusion matrix normalized to probabilities
%
% indivfold: Details of each individual fold. For mth run, kth fold,
% indivfold{dr,fe,cf} is a structure with these fields: 'confusion'
% is the confusion matrix of this fold. 'traintable' is a list of (1)
% indices and (2) their label indices (starts with 0) used for training
% this fold. 'testtable' is a table of four plus Nclass column vectors:
% (1) indices used for testing (correspond to the original sequence in
% Data and Labels, (2) truth label index (starts with 0),
% (3) decision label index (starts with 0) by classifier, (4) index
% (starts with 0) of the subspace that was used to reach the decision.
% (5...) Posterior probabilities
%
%
%
% USAGE EXAMPLE
%
% [pcorrect,pconfuse] = dataproc_main_multiintervalidation (Data,
% Labels, {'pca','cpca'}, {1,1}, {'aida','aida','aida'},
% {1,2,3}, {[],[],[],[]}, {'linear','quadratic','parzen','knn'},
% 'empirical', {0,0,0,3})
%
% Two dimension reduction algorithms (pca (parm 1) first, then cpca (parm
% 1)) are performed. Then three feature extraction are performed at
% different parameters (1,2,3) (they are all aida, but can be different if
% needed). Then four classifiers are performed using four corresponding
% parameters. The resulting output variables would be a cell in the format
% {DR,FE,CF}.
%
%
% AUTHOR
%
% Po T Wang
% 2009-10-08 (YYYY-MM-DD)
% 2009-10-09 (reduced memory usage and fixed LabelEnums bug)
% ptwang@uci.edu
% Biomedical Engineering, University of California Irvine
%


function [pcorrect, pconfuse, indivfold] = intervalidation (Data, Labels, TestIdx, DRfunC, DRparmC, FEfunC, FEparmC, CfunC, CnameC, Prior, CparmC, Opts)

Trialstotest = TestIdx;

%% Getting the parameters
Labels = Labels(:);
Nobs = length(Labels);
if Nobs ~= size(Data,1)
    error('Number of observations from Data and Labels disagree.');
end

TestIdx = unique(TestIdx(:).');

classes = unique(Labels);
Nclass = length(classes);
if Nclass <= 1
    error('At least two classes are required.');
end

% Convert labels to enumerated values
LabelEnums = Labels;
for c = 1:Nclass
    LabelEnums(find(Labels==classes(c)),1) = (c-1);
end

DRfunC = internalfunc_forcecell(DRfunC);
DRparmC = internalfunc_forcecell(DRparmC);
FEfunC = internalfunc_forcecell(FEfunC);
FEparmC = internalfunc_forcecell(FEparmC);
CfunC = internalfunc_forcecell(CfunC);
CnameC = internalfunc_forcecell(CnameC);
CparmC = internalfunc_forcecell(CparmC);

if length(DRfunC) ~= length(DRparmC)
    error('Lengths of DRfunC and DRparmC differ.');
elseif length(FEfunC) ~= length(FEparmC)
    error('Lengths of FEfunC and FEparmC differ.');
elseif length(CfunC) ~= length(CnameC) || length(CfunC) ~= length(CparmC)
    error('Lengths of CfunC, CnameC, and/or CparmC differ');
end

for dr = 1:length(DRfunC)
    if strcmp(DRfunC{dr},'cpca')
        DRfunC{dr} = @dataproc_func_cpca;
    elseif strcmp(DRfunC{dr},'pca')
        DRfunC{dr} = @dataproc_func_pca;
    elseif strcmp(DRfunC{dr},'ctpca')
        DRfunC{dr} = @dataproc_func_ctpca;
    end
end

for fe = 1:length(FEfunC)
    if strcmp(FEfunC{fe},'ida')
        FEfunC{fe} = @dataproc_func_ida;
    elseif strcmp(FEfunC{fe},'aida')
        FEfunC{fe} = @dataproc_func_aida;
    elseif strcmp(FEfunC{fe},'lda')
        FEfunC{fe} = @dataproc_func_lda;
    end
end

for cf = 1:length(CfunC)
    if length(CfunC{cf}) == 0
        CfunC{cf} = @dataproc_func_classify;
    end
end

if length(who('Opts')) == 0
    Opts = 1;
end

if nargout >= 3
    SWfindwinnersubspace = 1;
else
    SWfindwinnersubspace = 0;
end

if Opts
    SWprogressoutput = 1;
else
    SWprogressoutput = 0;
end

%% Partitioning data

for dr = 1:length(DRfunC)
    for fe = 1:length(FEfunC)
        for cf = 1:length(CfunC)
            confusion{dr,fe,cf} = zeros(Nclass,Nclass);
        end
    end
end

pcmax = 0;

% everything not test is train
TrainIdx = setdiff(1:size(Data,1),TestIdx);

for dr = 1:length(DRfunC)
    cvDatac = [];
    cvtestdatac = [];
    DRfun = DRfunC{dr};
    DRparm = DRparmC{dr};

    % Dimension reduction
    if length(DRfun)
        % DRmat can potentially be a cell. Always treat it as a cell
        DRmat = DRfun(Data(TrainIdx,:),LabelEnums(TrainIdx),DRparm);
        if ~iscell(DRmat)
            tmp = DRmat;
            clear('DRmat');
            DRmat{1} = tmp;
            clear('tmp');
        end
        for s = 1:length(DRmat)
            cvDatac{s} = Data(TrainIdx,:) * DRmat{s};
            cvtestdatac{s} = Data(TestIdx,:) * DRmat{s};
        end
    else
        DRmat = {};
        cvDatac{1} = Data(TrainIdx,:);
        cvtestdatac{1} = Data(TestIdx,:);
    end


    for fe = 1:length(FEfunC)

        cvfDatac = [];
        cvftestdatac = [];

        FEfun = FEfunC{fe};
        FEparm = FEparmC{fe};

        % Feature extraction
        if length(FEfun)
            FEmat{1} = FEfun(cvDatac{1},LabelEnums(TrainIdx),FEparm);
            cvfDatac{1} = cvDatac{1} * FEmat{1};
            cvftestdatac{1} = cvtestdatac{1} * FEmat{1};
            for s = 2:length(DRmat)
                FEmat{s} = FEfun(cvDatac{s},LabelEnums(TrainIdx),FEparm);
                cvfDatac{s} = cvDatac{s} * FEmat{s};
                cvftestdatac{s} = cvtestdatac{s} * FEmat{s};
            end
        else
            FEmat = {};
            cvfDatac = cvDatac;
            cvftestdatac = cvtestdatac;
        end


        for cf = 1:length(CfunC)
            Cfun = CfunC{cf};
            Cname = CnameC{cf};
            Cparm = CparmC{cf};

            % Classification

            if SWfindwinnersubspace
                [Decision,Posterior,WinnerSubspace] = Cfun(cvftestdatac,cvfDatac,LabelEnums(TrainIdx),Cname,Prior,Cparm);

            else
                Decision = Cfun(cvftestdatac,cvfDatac,LabelEnums(TrainIdx),Cname,Prior,Cparm);
            end
            Decision = Decision(:);

            for x = 0:Nclass-1
                for y = 0:Nclass-1
                    Tmat.confusion(x+1,y+1) = length(intersect(find(Decision == x), find(LabelEnums(TestIdx) == y)));
                    confusion{dr,fe,cf}(x+1,y+1) = confusion{dr,fe,cf}(x+1,y+1) + Tmat.confusion(x+1,y+1);
                end
            end

            if SWfindwinnersubspace
                Tmat.traintable = [TrainIdx(:) LabelEnums(TrainIdx)];
                Tmat.testtable = [TestIdx(:) LabelEnums(TestIdx) Decision WinnerSubspace-1 Posterior];
                indivfold{dr,fe,cf} = Tmat;
            end

        end % End of a classifier
    end % End of a feature extraction
end % End of a dimension reduction

for dr = 1:length(DRfunC)
    for fe = 1:length(FEfunC)
        for cf = 1:length(CfunC)
            pconfuse{dr,fe,cf}(:,:) = confusion{dr,fe,cf}(:,:) ./ (ones(size(confusion{dr,fe,cf}(:,:),1),1) * sum(confusion{dr,fe,cf}(:,:),1));
            Nactualclasstrial = sum(confusion{dr,fe,cf}(:,:),1);
            Ntotal = sum(Nactualclasstrial);

            if strcmp(Prior,'empirical')
                prior = Nactualclasstrial ./ sum(Nactualclasstrial);
            elseif strcmp(Prior,'equal')
                prior = ones(1,Nclass) ./ Nclass;
            else
                prior = Prior ./ sum(Prior);
            end
            pcorrect{dr,fe,cf} = sum(diag(pconfuse{dr,fe,cf}(:,:)).' .* prior);
        end % End of CF
    end % End of FE
end % End of DR
if SWprogressoutput
    for dr = 1:length(DRfunC)
        for fe = 1:length(FEfunC)
            for cf = 1:length(CfunC)
                if pcorrect{dr,fe,cf} > pcmax
                    pcmax = pcorrect{dr,fe,cf};
                end
            end
        end
    end
    fprintf(' Best Pcorrect = %-5.3f (chance = %-5.3f)\n',pcmax,1/Nclass);
end

function V = internalfunc_forcecell (V)
if ~iscell(V)
    tmp = V;
    clear('V');
    V{1} = tmp;
end

