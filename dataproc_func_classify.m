% [Decision, Posterior, WinningSubspaceID] = dataproc_func_classify (
% testdatac, traindatac, trainlabels, Classifier, Prior, Cparm )
%
% Classifies testdata using the classifier or pdf of your choice
%
% Decision: The class label that the classifier decides your test data
% belongs to.
%
% Posterior: The posterior probabilities among all classes.
%
% WinningSubspaceID: The index (starts from 1) of subspace, if applicable,
% that the decision is based upon. If only 1 subspace, returns 1.
%
% testdatac: Test data (each row is a different test data). The dimension
% (number of columns) should be fairly low, around 1~3. This parameter can
% be a cell containing test data sets from different subspaces. The number
% of subspaces must be identical to that of traindatac.
%
% traindatac: Training data (each row is a training sample). This parameter
% can be a cell containing training sets from different subspaces. However,
% the labels must still be identical. If it is a cell, highest posterior
% from each subspace is chosen and re-normalized among all subspaces.
%
% trainlabels: Class labels corresponding to the rows of training data.
%
% Classifier: Can be 'knn', 'wknn', 'parzen', 'quadratic', 'linear',
% 'diaglinear', 'diagquadratic'. Default: 'quadratic'. This is the
% probability density function or the likelihood function.
%
% Prior: Prior probabilities. Can be 'empirical', 'equal', or a horizontal
% vector. Default: 'empirical'
%
% Cparm: Classifier parameter. For knn and wknn, this is the number of
% neighbors. For parzen, this is the width multiplier. If Cparm == 0, it
% will be chosen automatically.
% Cparm: (advanced usage) For quadratic and linear, leave blank to use
% training data. If Cparm is a struct with the following structure,
% sufficient statistics mode is used:
%  Cparm.classes(c)   for the label name of cth class. This can be an
%                     integer or character
%  Cparm.stats{s,c,1} for the mean of sth subspace, cth class
%  Cparm.stats{s,c,2} for the variance of sth subspace, cth class
% If sufficient statistics mode is used, only quadratic or linear
% discriminants can be used, traindatac and trainlabels are ignored. If
% Prior other than 'equal' is needed, it must be specified numerically.
%
%
% Po T Wang
% 2008-08-22 (YYYY-MM-DD)
% 2009-05-12 - Sufficient stats
% ptwang@uci.edu
% University of California Irvine


function [Decision, Posterior, WinningSubspaceID] = dataproc_func_classify ( testdatac, traindatac, trainlabels, Classifier, Prior, Cparm )


if ~iscell(testdatac)
    tmp = testdatac;
    clear('testdatac');
    testdatac{1} = tmp;
end

if ~isempty(who('Cparm')) && isstruct(Cparm) && isfield(Cparm,'classes') && isfield(Cparm,'stats')
    SWsuffstat = 1;
else
    SWsuffstat = 0;
end

if SWsuffstat

    %20130522 somehow, clearing takes time %clear('traindatac','trainlabels');
    [Nsubspace,Nclass,tmp] = size(Cparm.stats); %#ok<NASGU>
    classes = Cparm.classes;
    if Nclass ~= length(Cparm.classes)
        error('Cparm.classes need to be the same length as number of classes');
    end
    Ndim = size(Cparm.stats{1,1,1}); %#ok<NASGU>

    % Assume Ndim is consistent for all subspaces and classes. (It should,
    % if the user has a clue, otherwise it won't work!)

    if isempty(who('Prior')) || isempty(Prior) || strcmpi(Prior,'equal')
        Prior = ones(1,Nclass) ./ Nclass;
    else
        Prior = Prior ./ sum(Prior);
    end
    
else

    if ~iscell(traindatac)
        tmp = traindatac;
        clear('traindatac');
        traindatac{1} = tmp;
    end

    Nsubspace = length(traindatac);
    Nobs = size(traindatac{1},1);
    Ndim = size(traindatac{1},2);

    for s = 2:Nsubspace
        if size(traindatac{s},2) ~= Ndim
            Ndim = min([Ndim,size(traindatac{s},2)]);
        end
    end

    if Nobs ~= length(trainlabels)
        error('Number of training samples from traindatac and trainlabels disagree.');
    end

    for s = 2:Nsubspace
        if size(traindatac{s},1) ~= Nobs
            error('One or more cell element in traindatac has different number of training samples.');
        end
        if size(traindatac{s},2) ~= Ndim
            %error('One or more cell element in traindatac has different dimension.');
            traindatac{s} = traindatac{s}(:,1:Ndim);
        end
    end

    classes = unique(trainlabels);
    Nclass = length(classes);

    if Nclass <= 1
        error('At least two unique labels are required.');
    end

    for c = 1:Nclass
        NtrainA(c) = length(find(trainlabels==classes(c))); %#ok<AGROW>
    end

    if isempty(who('Prior')) || strcmp(Prior,'empirical')
        Prior = NtrainA ./ Nobs;
    elseif isempty(Prior) || strcmp(Prior,'equal')
        Prior = ones(1,Nclass) ./ Nclass;
    end

end


if Nsubspace == 0
    error('At least one set of training is required.');
end

if isempty(who('Classifier'))
    Classifier = 'quadratic';
end

Posterior = nan(size(testdatac{1},1),Nclass,Nsubspace);
f = Posterior;
logf = Posterior;

switch Classifier
    case 'knn'
        if ~isempty(who('Cparm')) && ~isstruct(Cparm)
            Kn = Cparm;
        else
            Kn = [];
        end
%         for s = 1:Nsubspace
%             [Dec, Posterior(:,:,s)] = dataproc_func_knnclassify(testdatac{s},traindatac{s},trainlabels,Kn);
%         end
        
        for s = 1:Nsubspace
            N = 0;
            for c = 1:Nclass
                f(:,c,s) = dataproc_func_knnprob(testdatac{s},traindatac{s},trainlabels,Kn,classes(c));
                N = N + f(:,c,s) * Prior(c);

            end
            for c = 1:Nclass
                Posterior(:,c,s) = f(:,c,s) .* Prior(c) ./ N;
            end
        end
    case 'wknn'
        if ~isempty(who('Cparm')) && ~isstruct(Cparm)
            Kn = Cparm;
        else
            Kn = [];
        end
%         for s = 1:Nsubspace
%             [Dec, Posterior(:,:,s)] = dataproc_func_knnclassify(testdatac{s},traindatac{s},trainlabels,Kn);
%         end
        
        for s = 1:Nsubspace
            N = 0;
            for c = 1:Nclass
                f(:,c,s) = dataproc_func_knnprob(testdatac{s},traindatac{s},trainlabels,Kn,classes(c));
                N = N + f(:,c,s) * Prior(c);
            end
            for c = 1:Nclass
                Posterior(:,c,s) = f(:,c,s) .* Prior(c) ./ N;
            end
        end

    case 'parzen'
        if ~isempty(who('Cparm')) && ~isstruct(Cparm)
            Width = Cparm;
        else
            Width = [];
        end

        for s = 1:Nsubspace
            N = 0;
            for c = 1:Nclass
                %f(:,c,s) = dataproc_func_parzenwinpdf(testdatac{s},traindatac{s}(find(trainlabels==classes(c)),:),Width);
                f(:,c,s) = dataproc_func_parzenwinpdf(testdatac{s},traindatac{s}(trainlabels==classes(c),:),Width);
                N = N + f(:,c,s) * Prior(c);

            end
            for c = 1:Nclass
                Posterior(:,c,s) = f(:,c,s) .* Prior(c) ./ N;
            end
        end

    case 'quadratic'  % Our quadratic is 60% faster than MATLAB's
        for s = 1:Nsubspace
            N = 0;
            if SWsuffstat
                for c = 1:Nclass
                    f(:,c,s) = dataproc_func_normalpdf(testdatac{s},Cparm.stats{s,c,1},Cparm.stats{s,c,2});
                    N = N + f(:,c,s) * Prior(c);
                end
            else
                for c = 1:Nclass
                    %[f(:,c,s) logf(:,c,s)] = dataproc_func_normalpdf(testdatac{s},traindatac{s}(find(trainlabels==classes(c)),:));
                    [f(:,c,s) logf(:,c,s)] = dataproc_func_normalpdf(testdatac{s},traindatac{s}(trainlabels==classes(c),:));
                    N = N + f(:,c,s) * Prior(c);
                end
            end
            
            for c = 1:Nclass
                Posterior(:,c,s) = f(:,c,s) .* Prior(c) ./ N;
            end
            
            if any(N==0)
                [tmp I] = max(logf(:,:,s)+log(Prior(c)),[],2); %#ok<ASGLU>
                notI = [1:I-1,I+1:Nclass];
                Posterior(N==0,I,s) = 1;
                Posterior(N==0,notI,s) = 0;
            end
        end
    case 'linear'    % Our linear is 60% faster than MATLAB's
        for s = 1:Nsubspace
            N = 0;
            if SWsuffstat
                Sigma = Prior(1)*Cparm.stats{s,1,2};
                for c = 2:Nclass
                    Sigma = Sigma + Prior(c)*Cparm.stats{s,c,2};
                end
            else
                %Sigma = Prior(1)*cov(traindatac{s}(find(trainlabels==classes(1)),:));
                Sigma = Prior(1)*cov(traindatac{s}(trainlabels==classes(1),:));
                for c = 2:Nclass
                    %Sigma = Sigma + Prior(c)*cov(traindatac{s}(find(trainlabels==classes(c)),:));
                    Sigma = Sigma + Prior(c)*cov(traindatac{s}(trainlabels==classes(c),:));
                end
            end
            if SWsuffstat
                for c = 1:Nclass
                    f(:,c,s) = dataproc_func_normalpdf(testdatac{s},Cparm.stats{s,c,1},Sigma);
                    N = N + f(:,c,s) * Prior(c);
                end
            else
                for c = 1:Nclass
                    %[f(:,c,s), logf(:,c,s)] = dataproc_func_normalpdf(testdatac{s},traindatac{s}(find(trainlabels==classes(c)),:),Sigma);
                    [f(:,c,s), logf(:,c,s)] = dataproc_func_normalpdf(testdatac{s},traindatac{s}(trainlabels==classes(c),:),Sigma);
                    N = N + f(:,c,s) * Prior(c);
                end
            end
            for c = 1:Nclass
                Posterior(:,c,s) = f(:,c,s) .* Prior(c) ./ N;
            end
            if any(N==0)
                [tmp I] = max(logf(:,:,s)+log(Prior(c)),[],2); %#ok<ASGLU>
                notI = [1:I-1,I+1:Nclass];
                Posterior(N==0,I,s) = 1;
                Posterior(N==0,notI,s) = 0;
            end
        end
        
    otherwise
        for s = 1:Nsubspace
            [Dec, tmp, Posterior(:,:,s)] = classify(testdatac{s},traindatac{s},trainlabels,Classifier,Prior); %#ok<ASGLU>
        end

end

[Y,Isub] = max(Posterior,[],3);  % Main decider: Strongest posterior wins
Z = mean(Posterior,3); % Tie breaker if needed: Strongest mean wins
Posterior = Y ./ (sum(Y,2) * ones(1,Nclass));

% Detect ties (If winner is unique, I1 should == I2)
[Y,I1] = max(Posterior,[],2); %#ok<ASGLU>
[Y,I2] = max(fliplr(Posterior),[],2); %#ok<ASGLU>
I2 = Nclass+1-I2;

% Tie breaking
%Posterior(find(I1~=I2),:) = Z(find(I1~=I2),:);
Posterior(I1~=I2,:) = Z(I1~=I2,:);

% Decision
[Y,I] = max(Posterior,[],2); %#ok<ASGLU>
%WinningSubspaceID = diag(Isub(:,I)) % runs out of memory this way.
Isub = Isub.';
WinningSubspaceID = Isub((0:Nsubspace:Nsubspace*(size(Isub,2)-1))'+I);
Decision(:,1) = classes(I);
