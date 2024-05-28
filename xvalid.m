
% pcorrect and pconfuse are cells having dimensions {dr,fe,cf,pr}
% corresponding to DRFUN{dr}, FEFUN{fe}, CFUN{cf}, PRIOR{pr} that are
% specified by user or in the default xvalid settings

function [pcorrect, pconfuse, indivfold, ChanNames, TimePoints, FilterImage, DRFUN, FEFUN, FEPARM, CFUN, CNAME, PRIOR, BestParms, Fmat, DRPARM] = xvalid (TrainData, TrainLabels, varargin)

OPTS = 1;
KFOLD = 10;
MRUN = 10;
DRFUN = {'cpca'};
DRPARM = {[]};
FEFUN = {'aida','aida','aida','lda'};
FEPARM = {3,2,1,1};
CFUN = {'',''};
CNAME = {'quadratic','linear'};
PRIOR = {'empirical'};
CPARM = {[],[]};
SWsilent = 0;
sdfactor = 1;
COLUMN = 2;
COLOVERRIDE = 0;
ChanNames = {};
Nchan = [];
TimeRange = [];
TimeNames = {};
LabelNames = {};
thres = [];
FileName = [];
SWzscore = 0;
SWmzscore = 0;
SWlogscale = 0;
SWkeepempty = 1;
Identifier = '';
TRIALSTOTEST = [];
SWsortchans = 1;
BALANCED = 0;
BestParms = {};
NOPLOT = 0;
FPID = 1;
TimePoints = [];
FilterImage = {};
SWcompactfilter = 0;

if (~exist('TrainLabels', 'var') || isempty(TrainLabels)) && iscell(TrainData)
    Nclass = numel(TrainData);
    tmpTrainLabels = [];
    for c = 1:Nclass
        if size(TrainData{c},3) > 1
            n = size(TrainData{c},3);
            d = 3;
        else
            n = size(TrainData{c},1);
            d = 1;
        end
        tmpTrainLabels = [tmpTrainLabels; ones(n,1)*(c-1)];
    end
    TrainData = cat(d,TrainData{:});
    TrainLabels = tmpTrainLabels;
    clear tmpTrainLabels
end

TrainLabels = TrainLabels(:);


if ~isnumeric(TrainLabels) && iscell(TrainLabels) && length(TrainLabels) > 1 && ischar(TrainLabels{1})
    [TrainLabels, uniquelabels] = categorical_to_numeric_labels(TrainLabels);
    fprintf('The text labels have been converted to numerical labels:\n');
    for i = 1:length(uniquelabels)
        fprintf('\t%i\t%s\n', i, uniquelabels{i});
    end
end


if mod(length(varargin),2)
    error('Properties list must be in pairs, ie. property1 name, property1 value, ...');
end

for i = 1:2:length(varargin)
    switch upper(varargin{i})
        case 'DRFUN'
            if iscell(varargin{i+1}) && ~isempty(varargin{i+1}) && iscell(varargin{i+1}{1})
                % 20150316 Double-wrap fix
                DRFUN = lower(varargin{i+1}{1});
            else
                DRFUN = lower(varargin{i+1});
            end
        case 'DRPARM'
            if iscell(varargin{i+1}) && ~isempty(varargin{i+1}) && iscell(varargin{i+1}{1})
                % 20150316 Double-wrap fix
                DRPARM = varargin{i+1}{1};
            else
                DRPARM = varargin{i+1};
            end
        case 'FEFUN'
            if iscell(varargin{i+1}) && ~isempty(varargin{i+1}) && iscell(varargin{i+1}{1})
                % 20150316 Double-wrap fix
                FEFUN = lower(varargin{i+1}{1});
            else
                FEFUN = lower(varargin{i+1});
            end
            
        case 'FEPARM'
            if iscell(varargin{i+1}) && ~isempty(varargin{i+1}) && iscell(varargin{i+1}{1})
                % 20150316 Double-wrap fix
                FEPARM = varargin{i+1}{1};
            else
                FEPARM = varargin{i+1};
            end
        case 'KFOLD'
            KFOLD = varargin{i+1};
        case 'MRUN'
            MRUN = varargin{i+1};
        case 'PRIOR'
            PRIOR = varargin{i+1};
        case 'CNAME'
            CNAME = varargin{i+1};
        case 'DENSITY'
            CNAME = varargin{i+1};
        case 'CPARM'
            CPARM = varargin{i+1};
        case 'CFUN'
            CFUN = varargin{i+1};
        case 'TIME'
            TIMESPEC = varargin{i+1};
        case 'TIMESPEC'
            TIMESPEC = varargin{i+1};
        case 'CHAN'
            CHANSPEC = varargin{i+1};
        case 'CHANSPEC'
            CHANSPEC = varargin{i+1};
        case 'LABELRENAME'
            LABELRENAME = varargin{i+1};
        case 'LABELSPEC'
            LABELSPEC = varargin{i+1};
        case 'OPTS'
            OPTS = varargin{i+1};
        case 'QUICK'
            if varargin{i+1} == 1
                KFOLD = -5;
                MRUN = 1;
            end
        case 'LEAVEONEOUT'
            if varargin{i+1} == 1
                KFOLD = 0;
                MRUN = 1;
            end
        case 'INTRAINING'
            if varargin{i+1} == 1
                KFOLD = 1;
                MRUN = 1;
            end
        case 'SILENT'
            if varargin{i+1} == 1
                SWsilent = 1;
                OPTS = 0;
            end
        case 'COLUMN'
            if varargin{i+1} == 2
                COLUMN = 2;
                COLOVERRIDE = 1;
            elseif varargin{i+1} == 1
                COLUMN = 1;
                COLOVERRIDE = 1;
            end
        case 'COLUMNS'
            if varargin{i+1} == 2
                COLUMN = 2;
                COLOVERRIDE = 1;
            elseif varargin{i+1} == 1
                COLUMN = 1;
                COLOVERRIDE = 1;
            end
        case 'CHANNAMES'
            if ~isempty(varargin{i+1})
                ChanNames = varargin{i+1};
                if isempty(Nchan)
                    Nchan = length(ChanNames);
                end
            end
        case 'TIMENAMES'
            if ~isempty(varargin{i+1})
                TimeNames = varargin{i+1};
            end
        case 'TIMERANGE'
            if ~isempty(varargin{i+1})
                TimeRange = varargin{i+1};
            end
        case 'LABELNAMES'
            if ~isempty(varargin{i+1})
                LabelNames = varargin{i+1};
            end
        case 'NCHAN'
            if ~isempty(varargin{i+1})
                Nchan = varargin{i+1};
            end
        case 'FPRINTID'
            if ~isempty(varargin{i+1})
                FPID = varargin{i+1};
            end
        case 'THRES'
            if ~isempty(varargin{i+1})
                thres = varargin{i+1};
            end
        case 'FILENAME'
            if ~isempty(varargin{i+1})
                FileName = varargin{i+1};
            end
        case 'IDENTIFIER'
            if ~isempty(varargin{i+1})
                Identifier = varargin{i+1};
            end
        case 'ZSCORE'
            if varargin{i+1} == 1
                SWzscore = 1;
            elseif varargin{i+1} == 0
                SWzscore = 0;
            end
        case 'MZSCORE'
            if varargin{i+1} == 1
                SWmzscore = 1;
            elseif varargin{i+1} == 0
                SWmzscore = 0;
            end
        case 'LOGSCALE'
            if varargin{i+1} == 1
                SWlogscale = 1;
            elseif varargin{i+1} == 0
                SWlogscale = 0;
            end
        case 'KEEPEMPTY'
            if varargin{i+1} == 1
                SWkeepempty = 1;
            elseif varargin{i+1} == 0
                SWkeepempty = 0;
            end
        case 'TRIALSTOTEST'
            if ~isempty(varargin{i+1})
                TRIALSTOTEST = varargin{i+1};
            end
        case 'BALANCED'
            if varargin{i+1} == 2
                BALANCED = 2;
            elseif varargin{i+1} == 1
                BALANCED = 1;
            elseif varargin{i+1} == 0
                BALANCED = 0;
            end
        case 'NOPLOT'
            if varargin{i+1} == 1
                NOPLOT = 1;
            elseif varargin{i+1} == 0
                NOPLOT = 0;
            end
        case 'COMPACTIMAGE'
            if varargin{i+1} == 1
                SWcompactfilter = 1;
            elseif varargin{i+1} == 0
                SWcompactfilter = 0;
            end
    end
end

if length(DRFUN) > 1
    DRFUN = DRFUN(1);
    DRPARM = DRPARM(1);
    if OPTS > 0
        fprintf(FPID, 'Only the first dimension reduction function (DRFUN=%s) is used.\n', DRFUN{1});
    end
end
    


if (isempty(who('TrainLabels')) || isempty(TrainLabels)) && ~isempty(TrainData) && isstruct(TrainData(1))
    if OPTS > 0
        fprintf(FPID, 'Reconstructing rawdata... ');
    end
    [TrainData,TrainLabels,ChanNames,samplerate,datevector] = bsmspack_convert_to_rawdata(TrainData,'center');
    [S T N] = size(TrainData);
    if OPTS > 0
        fprintf(FPID, ' %i sensors, %i time points, %i trials.\n',S,T,N);
    end
else
    if size(TrainData,3) > 1
        [S T N] = size(TrainData);
    else
        [N T] = size(TrainData);
        S = 1;
    end
end

% if ~iscell(FEFUN) && iscell(FEPARM)
%     FEPARM = 1;
% end
% 
% if iscell(FEFUN) && ~iscell(FEPARM)
%     FEFUN = 'aida';
% end


Data = [];
Labels = [];

Ntime = size(TrainData,2);
if size(TrainData,3) == 1 && exist('Nchan', 'var') && ~isempty(Nchan) && isnumeric(Nchan)
    Ntime = Ntime / Nchan;
end

NeedToSpecTimeNames = 0;
if ~isempty(who('TIMESPEC')) || ~isempty(who('CHANSPEC'))
   if size(TrainData,3) == 1
       warning('TrainData is already reshaped. TIMESPEC is ignored, and CHANSPEC behaves differently.');
       TIMESPEC = [1:size(TrainData,2)];
       if ~isempty(who('CHANSPEC'))
           if islogical(CHANSPEC)
               CHANSPEC = find(CHANSPEC);
           end
           
           TrainData = TrainData(:,CHANSPEC);
           if ~isempty(ChanNames) && length(ChanNames) ~= length(CHANSPEC)
               ChanNames = ChanNames(CHANSPEC);
               Nchan = length(ChanNames);
           end
           
           if OPTS > 0
               disp('Subset of TrainData extracted using:');
               disp(['CHANSPEC = ' summarize_idx(CHANSPEC)]);
           end
       else
           CHANSPEC = [1:size(TrainData,1)];
       end
   else
       if isempty(who('TIMESPEC'))
           TIMESPEC = [1:size(TrainData,2)];
       end
       if isempty(who('CHANSPEC'))
           CHANSPEC = [1:size(TrainData,1)];
       end
       
       %20131101: Always convert to integer indexing
       if islogical(TIMESPEC)
           TIMESPEC = find(TIMESPEC);
       end
       if islogical(CHANSPEC)
           CHANSPEC = find(CHANSPEC);
       end
       
       TrainData = TrainData(CHANSPEC,TIMESPEC,:);
       if ~isempty(ChanNames) && length(ChanNames) ~= length(CHANSPEC)
           ChanNames = ChanNames(CHANSPEC);
       end
       
       NeedToSpecTimeNames = 1;
       
       if OPTS > 0
           disp('Subset of TrainData extracted using:');
           disp(['CHANSPEC = ' summarize_idx(CHANSPEC)]);
           disp(['TIMESPEC = ' summarize_idx(TIMESPEC)]);
       end
   end
end


if ~isempty(TimeRange)
    
    if exist('TIMESPEC','var')
        TIMESPEC2 = TIMESPEC;
    else
        TIMESPEC2 = size(TrainData,2);
    end
    if length(TimeRange) == length(TIMESPEC2)
        TimeNames = TimeRange;
    else
        TimeNames = linspace( TimeRange(1), TimeRange(end), Ntime+1 );
        TimeNames = TimeNames(1:end-1);
        TimeNames = string_to_cell(num2str(TimeNames),' ');
        % Keep only 20 points
        FullTimeNames = TimeNames;
        TimeNames(setdiff(1:length(TimeNames),1:length(TimeNames)/20:length(TimeNames))) = {''};
    end
end

if ~iscell(TimeNames)
    if ischar(TimeNames)
        TimeNames = string_to_cell(TimeNames, ', ');
    elseif isnumeric(TimeNames)
        %2021-11-03
        if size(TimeNames,2) == 2 && size(TimeNames,1) == Ntime
            tmp2 = cell(Ntime,1);
            for i = 1:size(TimeNames,1)
                tmp = sprintf('%.3g–%.3g', TimeNames(i,1), TimeNames(i,2));
                tmp2{i} = tmp;
            end
            TimeNames = tmp2;
            clear tmp tmp2
        else
            TimeNames = string_to_cell(num2str(TimeNames(:).'),' ');
        end
    end
end

if NeedToSpecTimeNames
    if ~isempty(TimeNames) && length(TimeNames) ~= length(TIMESPEC)
        TimeNames = TimeNames(TIMESPEC);
        if exist('FullTimeNames','var')
            FullTimeNames = FullTimeNames(TIMESPEC);
        end
    end
end


LabelIDs = unique(TrainLabels(:).');
if isempty(LabelNames)
    LabelNames = string_to_cell(num2str(LabelIDs), ' ');
end
LabelIDtoNameMap(:,2) = LabelNames(:);
for i = 1:length(LabelIDs)
    LabelIDtoNameMap{i,1} = LabelIDs(i);
end


if SWlogscale
    
    if any(TrainData(:) <= 0)
        TrainData = eps+abs(TrainData);
        warning('TrainData contains non-positive numbers. To compute log, it has been fixed with abs and eps.');
    end
    TrainData = log(TrainData);
    if OPTS > 0
        disp('TrainData time/freq series has been converted to natural-log scale.');
    end
end

if SWzscore || SWmzscore
    if size(TrainData,3) == 1
        warning('TrainData is already reshaped. ZSCORE/MZSCORE directive is ignored!');
        SWzscore = 0;
        SWmzscore = 0;
    else
        if SWmzscore && SWzscore
            warning('Both ZSCORE and RZSCORE are enabled. Disabling ZSCORE and keeping MZSCORE.');
            SWzscore = 0;
        end
        
        if SWmzscore
            try
                TrainData = mzscore(TrainData,[],2);
            catch
                for ch = 1:size(TrainData,1)
                    for tr = 1:size(TrainData,3)
                        TrainData(ch,:,tr) = mzscore(TrainData(ch,:,tr));
                    end
                end
            end
        elseif SWzscore
            try
                TrainData = zscore(TrainData,[],2);
            catch
                for ch = 1:size(TrainData,1)
                    for tr = 1:size(TrainData,3)
                        TrainData(ch,:,tr) = zscore(TrainData(ch,:,tr));
                    end
                end
            end
        end
        if OPTS > 0
            disp('TrainData time/freq series has been normalized to z-scores.');
        end
    end
end

if ~isempty(who('LABELRENAME')) && length(LABELRENAME) >= 2
    % Correct format of LABELRENAME is pairs of old-new labels
    % [1old 1new 2old 2new ...]
    LabelRenameHistory = zeros(floor(length(LABELRENAME)/2),2);
    for i = 1:2:length(LABELRENAME)-1
        tmp = find(TrainLabels==LABELRENAME(i));
        tmp2 = find(TrainLabels==LABELRENAME(i+1));
        LabelRenameHistory((i+1)/2,:) = [length(tmp) length(tmp2)];
        TrainLabels(tmp) = LABELRENAME(i+1);
        if OPTS > 0
            disp(['Renamed label "' num2str(LABELRENAME(i)) '" (' num2str(length(tmp)) ' trials) to "' num2str(LABELRENAME(i+1)) '"']);
        end
    end
else
    LABELRENAME = [];
    LabelRenameHistory = [];
end

classes = unique(TrainLabels);
for j = 1:length(classes)
    tmp = find(TrainLabels==classes(j));
    if length(tmp) < 2
        keepidx = setdiff(1:length(TrainLabels(:)),tmp);
        TrainLabels = TrainLabels(keepidx);
        if size(TrainData,3) == 1
            TrainData = TrainData(keepidx,:);
        else
            TrainData = TrainData(:,:,keepidx);
        end
        if OPTS > 0
            disp(['Deleted class ' num2str(classes(j)) ' (only 1 trial)']);
        end
    end
end

if ~isempty(who('LABELSPEC'))
    LABELSPEC = unique(LABELSPEC(:)).';
    for i = 1:length(LABELSPEC)
        tridx{i} = find(TrainLabels==LABELSPEC(i));
    end
    atridx = cat(1,tridx{:});
    TrainLabels = TrainLabels(atridx);
    if size(TrainData,3) == 1
        TrainData = TrainData(atridx,:);
    else
        TrainData = TrainData(:,:,atridx);
    end
    if OPTS > 0
        disp('Only the following classes are included:');
        for i = 1:length(LABELSPEC)
            disp(['LABEL = ' num2str(LABELSPEC(i)) ' (' num2str(length(tridx{i})) ' trials)']);
        end
    end
end


if size(TrainData,3) > 1
    if OPTS > 0
        %disp('Reshaping data...');
    end
    Nchan = size(TrainData,1);
    TrainData = rdreshape(TrainData);
end

% if size(TrainLabels,2) > 1
%     TrainLabels = TrainLabels.';
% end

% if ~iscell(FEPARM) && FEPARM >= size(TrainData,2)
%     FEFUN = '';
%     FEPARM = [];
% end

Data = TrainData;
Labels = TrainLabels;

% % Add noise (for digital data)
% Data = Data ;%+ randn(size(Data))*(max(max(Data))-min(min(Data)))*1e-9;

%[pcorrect,pconfuse] = dataproc_main_crossvalidation(Data,Labels,KFOLD,MRUN,DRFUN,DRPARM,FEFUN,FEPARM,CFUN,CNAME,PRIOR,CPARM,OPTS);


classes = unique(Labels);
Nclass = length(classes);

% 20150316: Automatically add more dimensions for higher number of classes
if Nclass > 2 && nargin >= 4 && min(cellfun(@isempty,strfind(varargin(1:2:end), 'FEFUN'))) && min(cellfun(@isempty,strfind(varargin(1:2:end), 'FEPARM')))
    MAX_AIDA_DIM = max(3,Nclass - 1);
    FEFUN_AIDA = cell(1, MAX_AIDA_DIM);
    FEPARM_AIDA = cell(1, MAX_AIDA_DIM);
    for i = 1:MAX_AIDA_DIM
        FEFUN_AIDA{i} = 'aida';
        FEPARM_AIDA{i} = MAX_AIDA_DIM+1-i;
    end
    MAX_LDA_DIM = Nclass - 1;
    FEFUN_LDA = cell(1, MAX_LDA_DIM);
    FEPARM_LDA = cell(1, MAX_LDA_DIM);
    for i = 1:MAX_LDA_DIM
        FEFUN_LDA{i} = 'lda';
        FEPARM_LDA{i} = MAX_LDA_DIM+1-i;
    end
    FEFUN = [FEFUN_AIDA FEFUN_LDA];
    FEPARM = [FEPARM_AIDA FEPARM_LDA];
end

% 20120828: Always request indivfold. We will do fun stuff with it in xvalid.
[pcorrect, pconfuse, indivfold, Parameters] = dataproc_main_multicrossvalidation(Data,Labels,KFOLD,MRUN,DRFUN,DRPARM,FEFUN,FEPARM,CFUN,CNAME,PRIOR,CPARM,OPTS,TRIALSTOTEST,BALANCED);

if ~COLOVERRIDE
    if Nclass > 3
        COLUMN = 1;
    end
end

KFOLD = Parameters.K;
MRUN = Parameters.M;
NtrialA = Parameters.NtrialA;
TestCandidates = Parameters.TestCandidates;
SWdisablestratify = Parameters.SWdisablestratify;
PriorC = Parameters.PriorC;
BALANCED = Parameters.SWbalanced;
for i = 1:length(PriorC)
    switch PriorC{i}
        case 'equal'
            PriorC{i} = ones(1,Nclass)./Nclass;
        case 'empirical'
            PriorC{i} = NtrialA ./ sum(NtrialA);
    end
end




if MRUN == 1
    % If only one run of CV was performed, the std reported will be
    % between-folds (for k-fold) or between trials (for leave-one-out)
    %testtable = cell(size(indivfold));
    confusion = cell(size(indivfold));
    pcorrect = confusion;
    for dr = 1:size(indivfold,1)
        for fe = 1:size(indivfold,2)
            for cf = 1:size(indivfold,3)
                for pr = 1:size(indivfold,4)
                    m = 1;
                    tmp = cell(1,KFOLD);
                    tmp2 = tmp;
                    %tmp3 = zeros(KFOLD,6);
                    for k = 1:size(indivfold{dr,fe,cf,pr},2)
                        %tmp3(k,:) = indivfold{dr,fe,cf,pr}{m,k}.testtable;
                        tmp{k} = indivfold{dr,fe,cf,pr}{m,k}.confusion;
                        tmp2{k} = sum(diag(tmp{k})) / sum(tmp{k}(:));
                    end
                    %testtable{dr,fe,cf,pr} = tmp3;
                    confusion{dr,fe,cf,pr} = cat(3,tmp{:});
                    pcorrect{dr,fe,cf,pr} = cat(2,tmp2{:});
                    
                    
                    tmp_mu = mean(confusion{dr,fe,cf,pr},3);
                    tmp_fac = ones(Nclass,1)*sum(tmp_mu,1);
                    %tmp_sd = std(confusion{dr,fe,cf,pr},[],3);
                    %tmp_mu = tmp_mu ./ tmp_fac;
                    %tmp_sd = tmp_sd ./ tmp_fac;
                    
                    
                    
                    
                    for k = 1:size(confusion{dr,fe,cf,pr},3)
                        tmp = confusion{dr,fe,cf,pr}(:,:,k);
                        pconfuse{dr,fe,cf,pr}(:,:,k) = tmp ./ tmp_fac;
                        
                    end
                end
            end
        end
    end
end






for j = 1:Nclass
    Ntrial(j) = length(find(Labels==classes(j)));
end

if ~iscell(PRIOR)
    PRIOR = {PRIOR};
end
PR = length(PRIOR);
for pr = 1:PR
    Prior = PRIOR{pr};
    if ~SWsilent
        fprintf(FPID, '\n');
        %     if strcmpi(DRFUN,'cpca')
        %         cpcamat = dataproc_func_cpca(Data,Labels);
        %         for i = 1:length(cpcamat)
        %             cpcadim(i) = size(cpcamat{i},2);
        %         end
        %     elseif strcmpi(DRFUN,'pca')
        %         pcamat = dataproc_func_pca(Data,Labels);
        %         pcadim = size(pcamat,2);
        %     end
        
        if SWdisablestratify
            tmp = 'Unstratified';
        else
            tmp = 'Stratified';
        end
        if KFOLD == length(TestCandidates) || KFOLD == 0
            % This is leave-one-out
            fprintf(FPID, 'Leave-one-out Cross Validation Report\n');
            CANHAVEUNCERTAINTY = 2;
        elseif abs(KFOLD) == 1
            fprintf(FPID, '%s In-training Self Validation Report\n',tmp);
            CANHAVEUNCERTAINTY = 0;
        elseif MRUN == 1
            fprintf(FPID, '%s %i-fold Cross Validation Report\n',tmp,KFOLD);
            CANHAVEUNCERTAINTY = 2;
        else
            fprintf(FPID, '%s %i-fold Cross Validation Report\n',tmp,KFOLD);
            CANHAVEUNCERTAINTY = 1;
        end
        
        if ~isempty(LABELRENAME)
            fprintf(FPID, 'The following classes have been renamed or merged:\n')
            for i = 1:2:length(LABELRENAME)-1
                if LabelRenameHistory((i+1)/2,2) > 0
                    RenameType = ['merged with the existing ' num2str(LabelRenameHistory((i+1)/2,2)) ' trials'];
                else
                    RenameType = 'renamed';
                end
                fprintf(FPID, ' Class %i (%i trials) -> Class %i (%s)\n', LABELRENAME(i), LabelRenameHistory((i+1)/2,1), LABELRENAME(i+1), RenameType);
            end
            fprintf(FPID, 'The remainder of this report will use the new class names.\n');
        end
        fprintf(FPID, '  Number of total trials: ');
        for i = 1:Nclass-1
            fprintf(FPID, '%i (class %i), ',Ntrial(i),classes(i));
        end
        fprintf(FPID, '%i (class %i), %i (total)\n',Ntrial(Nclass),classes(Nclass), sum(Ntrial));
        fprintf(FPID, ' Number of tested trials: ');
        for i = 1:Nclass-1
            fprintf(FPID, '%i (class %i), ',Parameters.NtestcandidateA(i),classes(i));
        end
        fprintf(FPID, '%i (class %i), %i (total)\n',Parameters.NtestcandidateA(Nclass),classes(Nclass), sum(Parameters.NtestcandidateA));
        if BALANCED >= 2
            fprintf(FPID, 'BALANCED: Number of training and testing trials were balanced during each run.\n');
        elseif BALANCED >= 1
            fprintf(FPID, 'BALANCED: Number of training trials were balanced during each run.\n');
        end
        fprintf(FPID, '     Data dimension (starting): %i\n', size(Data,2));
        
        switch DRFUN{1}
            case 'pca'
                tmp = dataproc_func_pca(TrainData, TrainLabels, DRPARM{1});
                if iscell(tmp)
                    for tmp2 = 1:length(tmp)
                        drdim(tmp2) = size(tmp{tmp2},2);
                    end
                else
                    drdim = size(tmp,2);
                end
            case 'cpca'
                tmp = dataproc_func_cpca(TrainData, TrainLabels, DRPARM{1});
                if iscell(tmp)
                    for tmp2 = 1:length(tmp)
                        drdim(tmp2) = size(tmp{tmp2},2);
                    end
                else
                    drdim = size(tmp,2);
                end
            case 'fkt'
                tmp = dataproc_func_fktdr(TrainData, TrainLabels, DRPARM{1});
                if iscell(tmp)
                    for tmp2 = 1:length(tmp)
                        drdim(tmp2) = size(tmp{tmp2},2);
                    end
                else
                    drdim = size(tmp,2);
                end
            otherwise
                drdim = NaN;
                
        end
        
        fprintf(FPID, '     Data dimension (after DR): %s (estimated)\n', num2str(drdim));
        if ~isempty(ChanNames)
            fprintf(FPID, 'Included %i channels: %s\n', length(ChanNames), cell_to_string(ChanNames(:).', ', '));
        end
        if ~isempty(TimeNames)
            if exist('FullTimeNames','var') && ~isempty(FullTimeNames) && iscell(FullTimeNames)
                fprintf(FPID, 'Included %i time/freq: %s\n', length(TimeNames), cell_to_string(FullTimeNames(:).', ', '));
            else
                fprintf(FPID, 'Included %i time/freq: %s\n', length(TimeNames), cell_to_string(TimeNames(:).', ', '));
            end
        end
        if exist('LABELSPEC','var')
            [~, ~, ib] = intersect(LABELSPEC, LabelIDs);
            fprintf(FPID, 'Included %i classes: %s\n', length(LABELSPEC), cell_to_string(LabelNames(ib), ', '));
            %fprintf(FPID, 'Included %i classes: %s\n', length(LABELSPEC), regexprep(num2str(LABELSPEC(:).'),' +', ', '));
        end
        fprintf(FPID, 'Included %i test trials: %s\n', length(TestCandidates), get_contig_groups_string(TestCandidates));
        if ~iscell(DRPARM)
            DRPARM = {DRPARM};
        end
        fprintf(FPID, 'Dimension reduction: %s (%g)\n',DRFUN{1},DRPARM{1});
        %     if strcmpi(DRFUN,'cpca')
        %         fprintf(FPID, ' Reduced dimensions: ');
        %         for i = 1:length(cpcadim)-1
        %             fprintf(FPID, '%i (class %i), ',cpcadim(i),classes(i));
        %         end
        %         fprintf(FPID, '%i (class %i)\n',cpcadim(end),classes(end));
        %     elseif strcmpi(DRFUN,'pca')
        %         fprintf(FPID, '  Reduced dimension: ');
        %         fprintf(FPID, '%i\n',pcadim);
        %     end
        fprintf(FPID, 'Prior probabilities: [%s]\n',num2str(Prior));
        %disp([num2str(KFOLD) '-fold Cross Validation Report']);
        %disp(['       Class labels: ' num2str(classes)]);
        %disp(['# samples per class: ' num2str(Ntrial)]);
        %disp(['     Data dimension: ' num2str(size(Data,2))]);
        %disp(['Dimension reduction: ' DRFUN '(' num2str(DRPARM) ')']);
        %disp(['    CPCA dimensions: ' num2str(cpcadim)]);
        %disp(['   Density function: ' CNAME '(' num2str(CPARM) ')']);
        %disp(['Prior probabilities: ' '[' num2str(PRIOR) ']']);
        if CANHAVEUNCERTAINTY == 1
            fprintf(FPID, ' K-fold Uncertainty: ±%g stdev between the %i runs.\n',sdfactor, MRUN);
        elseif CANHAVEUNCERTAINTY == 2
            fprintf(FPID, ' K-fold Uncertainty: ±%g stdev between the %i folds.\n',sdfactor, KFOLD);
        elseif abs(KFOLD) == 1
            fprintf(FPID, '        Overfitting: Likely\n');
            %disp(['        Overfitting: Likely']);
        elseif KFOLD == length(TestCandidates) || KFOLD == 0
            % nothing here
        elseif MRUN == 1
            fprintf(FPID, ' Only one run. Stdev not available.\n');
        end
    end
    
    % % old system
    % pcormean = mean(cat(1,pcorrect{:}),2);
    % pcorstd = std(cat(1,pcorrect{:}),[],2);
    % pconmean = mean(cat(4,pconfuse{:}),3);
    % pconstd = std(cat(4,pconfuse{:}),[],3);
    
    % We do not test DRFUN in xvalid. Always use the first one.
    pcor = shiftdim(pcorrect(1,:,:,pr),1);
    pcon = shiftdim(pconfuse(1,:,:,pr),1);
    
    % i is for FE dimension
    % j=1 for quadratic, j=2 for linear
    for i = 1:size(pcor,1)
        for j = 1:size(pcor,2)
            pcormean(i,j) = mean(pcor{i,j});
            pcorstd(i,j) = std(pcor{i,j});
            pconmean(:,:,i,j) = mean(pcon{i,j},3);
            pconstd(:,:,i,j) = std(pcon{i,j},[],3);
        end
    end
    
    
    pcormean2(:,:,pr) = pcormean;
    
    
    %[Y,I] = sort(max(pcormean,[],2));
    %clear('Y');
    %I = flipud(I(:));
    
    [Y,I] = sort(reshape(pcormean,[],1),1,'descend');
    [Y J] = sort(I,1,'ascend');
    if numel(J) >= 2 && mod(numel(J),2) == 0 && length(CNAME) == 2
        J = reshape(J,[],2);
    end
    
    if ~SWsilent
        if COLUMN == 2
            fprintf(FPID, '-----------------------------------');
            for i = 3:Nclass
                fprintf(FPID, '---------------');
            end
            fprintf(FPID, '+------------------------------------');
            for i = 3:Nclass
                fprintf(FPID, '---------------');
            end
            fprintf(FPID, '+\n');
        end
        for fi = 1:length(FEPARM)
            if COLUMN == 2
                for cf = 1:length(CNAME)
                    RANKTEXT = ['(#' num2str(J(fi,cf)) ')'];
                    fprintf(FPID, '%5.5s(%i) %-4.4s, ',FEFUN{fi},FEPARM{fi},CNAME{cf});
                    if ~CANHAVEUNCERTAINTY
                        fprintf(FPID, '%6.1f%%      %6s',pcormean(fi,cf)*100,RANKTEXT);
                    else
                        fprintf(FPID, '%6.1f±%4.1f%% %6s',pcormean(fi,cf)*100,pcorstd(fi,cf)*100*sdfactor,RANKTEXT);
                    end
                    for i = 3:Nclass
                        fprintf(FPID, '               ');
                    end
                    fprintf(FPID, ' | ');
                    %disp('Confusion matrix');
                end
                
                fprintf(FPID, '\n');
                
                if ~CANHAVEUNCERTAINTY
                    onelen = 15;
                else
                    onelen = 15;
                end
                for cf = 1:length(CNAME)
                    if BALANCED
                        tmp = infotransferrate(eye(length(PriorC{pr})), ones(1,Nclass)/Nclass);
                    else
                        tmp = infotransferrate(eye(length(PriorC{pr})), PriorC{pr});
                    end
                    fprintf(FPID, ['%-' num2str(4+onelen*Nclass) 's'], sprintf('ITR= %.3g b/tr (%.3g max)', infotransferrate(pconmean(:,:,fi,cf), PriorC{pr}), tmp));
                    fprintf(FPID, ' | ');
                end
                fprintf(FPID, '\n');
                
                me = pconmean(:,:,fi,1)*100;
                sd = pconstd(:,:,fi,1)*100;
                for i = 1:size(me,1)
                    for cf = 1:length(CNAME)
                        me = pconmean(:,:,fi,cf)*100;
                        sd = pconstd(:,:,fi,cf)*100;
                        for j = 1:size(me,2)
                            if ~CANHAVEUNCERTAINTY
                                fprintf(FPID, '%6.1f%%        ', me(i,j));
                            else
                                fprintf(FPID, '%6.1f±%5.1f%%  ', me(i,j), sdfactor*sd(i,j));
                            end
                        end
                        fprintf(FPID, '     | ');
                    end
                    fprintf(FPID, '\n');
                end
                
                
                fprintf(FPID, '-----------------------------------');
                for i = 3:Nclass
                    fprintf(FPID, '---------------');
                end
                fprintf(FPID, '+------------------------------------');
                for i = 3:Nclass
                    fprintf(FPID, '---------------');
                end
                fprintf(FPID, '+\n');
            else
                for cf = 1:length(CNAME)
                    RANKTEXT = ['(Rank ' num2str(J(fi,cf)) ')'];
                    fprintf(FPID, '%s(%i), %s, ',FEFUN{fi},FEPARM{fi},CNAME{cf});
                    if ~CANHAVEUNCERTAINTY
                        fprintf(FPID, 'Pcorrect =%6.1f%% %s\n',pcormean(fi,cf)*100,RANKTEXT);
                        %disp([' FE: ' FEFUN{fi} '(' num2str(FEPARM{fi}) ')' ' ' 'CF: ' CNAME{cf} ': Pcorrect = ' num2str(pcormean(fi,cf)*100,'%5.1f') '%' ' ' RANKTEXT]);
                    else
                        fprintf(FPID, 'Pcorrect =%6.1f±%.1f%% %s\n',pcormean(fi,cf)*100,pcorstd(fi,cf)*100*sdfactor,RANKTEXT);
                        %fprintf(FPID, 'FE: %s(%i)\tCF: %s\tPcorrect = %5.1f%% (sd=%5.1f%%)\t%s\n',FEFUN{fi},FEPARM{fi},CNAME{cf},pcormean(fi,cf)*100,pcorstd(fi,cf)*100,RANKTEXT);
                        %disp([' Feature extraction: ' FEFUN{fi} '(' num2str(FEPARM{fi}) ')' ' ' 'Discriminant: ' CNAME{cf} ': Pcorrect = ' num2str(pcormean(fi,1)*100,'%5.1f') '% (sd=' num2str(pcorstd(fi,cf)*100,'%5.1f') '%)' ' ' RANKTEXT]);
                    end
                    if BALANCED
                        tmp = infotransferrate(eye(length(PriorC{pr})), ones(1,Nclass)/Nclass);
                    else
                        tmp = infotransferrate(eye(length(PriorC{pr})), PriorC{pr});
                    end
                    fprintf(FPID, '%-34s\n', sprintf('ITR= %.3g bits/trial (%.3g max)', infotransferrate(pconmean(:,:,fi,cf), PriorC{pr}), tmp));
                    %disp('Confusion matrix');
                    me = pconmean(:,:,fi,cf)*100;
                    sd = pconstd(:,:,fi,cf)*100;
                    if ~CANHAVEUNCERTAINTY
                        fprintf(FPID, '%8s', 'Dec\Tru ');
                    else
                        fprintf(FPID, '%13s', 'Dec\Tru ');
                    end
                    for j = 1:size(me,2)
                        if ~CANHAVEUNCERTAINTY
                            fprintf(FPID, '%8s', sprintf('[%i]', classes(j)));
                        else
                            fprintf(FPID, '%13s', sprintf('[%i]', classes(j)));
                        end
                    end
                    fprintf(FPID, '\n');
                    for i = 1:size(me,1)
                    if ~CANHAVEUNCERTAINTY
                        fprintf(FPID, '%8s', sprintf('[%i]', classes(i)));
                    else
                        fprintf(FPID, '%13s', sprintf('[%i]', classes(i)));
                    end
                        for j = 1:size(me,2)
                            if ~CANHAVEUNCERTAINTY
                                fprintf(FPID, '%8s', sprintf('%6.1f%%', me(i,j)));
                            else
                                fprintf(FPID, '%13s', sprintf('%6.1f±%5.1f%%', me(i,j), sdfactor*sd(i,j)));
                            end
                        end
                        fprintf(FPID, '\n');
                    end
                    fprintf(FPID, '\n');
                end % End CF
            end % End if COLUMN == 2
        end % End FEPARM
        fprintf(FPID, '\n');
    end % End if silent
end % End prior

% Find the bestest (restricts to CPCA and equal prior)
mm = max(max(pcormean2(:,:,1)));
for fe = 1:length(FEPARM)
    for cf = 1:length(CFUN)
        for pr = 1
            if pcormean2(fe,cf,pr) == mm
                festar = fe;
                cfstar = cf;
                prstar = pr;
                BestParms.DRFUN = DRFUN{1};
                BestParms.DRPARM = DRPARM{1};
                BestParms.FEFUN = FEFUN{festar};
                BestParms.FEPARM = FEPARM{festar};
                BestParms.CFUN = CFUN{cfstar};
                BestParms.CNAME = CNAME{cfstar};
                BestParms.CPARM = CPARM{cfstar};
                BestParms.PRIOR = PRIOR{1};
                BestParms.PCORRECT = pcorrect{1,fe,cf,pr};
                BestParms.PCONFUSE = pconfuse{1,fe,cf,pr};
                
                switch lower(BestParms.DRFUN)
                    case 'pca'
                        DRF = @dataproc_func_pca;
                    case 'cpca'
                        DRF = @dataproc_func_cpca;
                    case 'fkt'
                        DRF = @dataproc_func_fktdr;
                    otherwise
                        DRF = BestParms.DRFUN;
                end
                
                switch lower(BestParms.FEFUN)
                    case 'aida'
                        FEF = @dataproc_func_aida;
                    case 'lda'
                        FEF = @dataproc_func_lda;
                    otherwise
                        FEF = BestParms.FEFUN;
                end
                
                tm1 = DRF(Data,Labels,BestParms.DRPARM);
                Ftrain = {};
                tm2 = {};
                trfmatbest = {};
                if iscell(tm1)
                    for i = 1:length(tm1)
                        tm2{i} = FEF(Data*tm1{i},Labels,BestParms.FEPARM);
                        trfmatbest{i} = tm1{i} * tm2{i};
                        Ftrain{i} = Data*trfmatbest{i};
                    end
                else
                    tm2 = FEF(Data*tm1,Labels,BestParms.FEPARM);
                    trfmatbest = tm1 * tm2;
                    Ftrain{1} = Data*trfmatbest;
                    tm1 = {tm1};
                end
                
                BestParms.trfmat = trfmatbest;
                SuffStat.classes = unique(Labels);
                for s = 1:length(tm1)
                    for c = 1:length(SuffStat.classes)
                        SuffStat.stats{s,c,1} = mean(Ftrain{s}(Labels==SuffStat.classes(c),:),1);
                        SuffStat.stats{s,c,2} = cov(Ftrain{s}(Labels==SuffStat.classes(c),:));
                    end
                end
                BestParms.SuffStat = SuffStat;
                
                break
            end
        end
    end
end


if ~SWsilent
    if Nclass > 3
        fprintf(FPID, 'Result tables are very large. Scroll up to see both.\n');
    end
end

% Generate filter image now if requested

if ~isempty(Nchan)
    % Determine which 1D classifier has highest performance
    ww = -1;
    ii = -1;
    for fe = 1:length(FEFUN)
        for cf = 1:length(CNAME)
            score = pcormean(fe,cf)-pcorstd(fe,cf)/2;
            if FEPARM{fe} == 1 && score > ww
                ww = score;
                ii = fe;
            end
        end
    end
    if ii == -1
        % No 1D classifier was specified
        for fe = 1:length(FEFUN)
            for cf = 1:length(CNAME)
                score = pcormean(fe,cf)-pcorstd(fe,cf)/2;
                if score > ww
                    ww = score;
                    ii = fe;
                end
            end
        end
    end
    switch FEFUN{ii}
        case 'aida'
            trfmat = cpca_aida(Data,Labels,1);
            ClassifierDesc = 'CPCA+AIDA';
        case 'lda'
            trfmat = cpca_lda(Data,Labels,1);
            ClassifierDesc = 'CPCA+LDA';
        case 'ida'
            trfmat = cpca_ida(Data,Labels,1);
            ClassifierDesc = 'CPCA+IDA';
        otherwise
            trfmat = [];
    end
    for s = 1:length(trfmat)
        Fmat{s} = reshape(trfmat{s}(:,1),Nchan,[]);
        %Classname{s} = num2str(classes(s));
        Classname{s} = LabelNames{classes(s) == LabelIDs};
    end
    impact = zeros(length(Fmat),Nchan);
    for s = 1:length(Fmat)
        for i = 1:Nchan
            impact(s,i) = Fmat{s}(i,:)*Fmat{s}(i,:).'/Nchan^2;
        end
    end
    impact = mean(impact,1);
    
    if ~isempty(ChanNames)
        FmatChanNames = ChanNames;
        FmatNchan = Nchan;
        %trfmat_trimmed = trfmat;
        if SWkeepempty ~= 1
            ToDelete = find(impact<1e-7);
            trfmat_trimmed = cell(1,length(Fmat));
            for s = 1:length(Fmat)
                Fmat{s} = Fmat{s}(setdiff(1:Nchan,ToDelete),:);
                %trfmat_trimmed{s} = reshape(Fmat{s},[],1);
            end
            FmatChanNames = FmatChanNames(setdiff(1:Nchan,ToDelete));
            FmatNchan = length(FmatChanNames);
            if OPTS > 0 && ~isempty(ToDelete)
                disp('Channels with no class differences will not be displayed in the heat map.');
                disp(['These channels are therefore excluded: ' cell_to_string(ChanNames(ToDelete),', ')]);
            end
        end
        
        
        clear tmp_pcor tmp_pcondiag CommentTextCell
        tmp_pcor = [num2str(pcormean(J == min(J(:)))*100,'%.1f') 177 num2str(pcorstd(J == min(J(:)))*100,'%.1f%%')];
        tmp_Ntrial = Ntrial;
        if BALANCED
            tmp_Ntrial(:) = min(tmp_Ntrial);
        end
        for c = 1:size(pconmean,1)
            tmp_pcondiag{c} = [num2str(pconmean(c,c,J == min(J(:)))*100,'%.1f') 177 num2str(pconstd(c,c,J == min(J(:)))*100,'%.1f%%')];
            if SWzscore
                tmp_CmtZscored = 'Z-scored, ';
            else
                tmp_CmtZscored = '';
            end
            if isempty(Identifier)
                Identifier = 'NaN';
            end
            CommentTextCell(c) = {['ID=' Identifier ', ' tmp_CmtZscored 'P[all]=' tmp_pcor ', P[' Classname{c} ']=' tmp_pcondiag{c} ', n_' num2str(c-1) '=' num2str(tmp_Ntrial(c)) ', nch=' num2str(FmatNchan)]};
        end
        
        if SWcompactfilter
            if ~isempty(FileName)
                generate_compact_filter(trfmat, FmatChanNames, [], TimeNames, FileName);
                FilterImage = [];
                TimePoints = [];
            else
                generate_compact_filter(trfmat, FmatChanNames, [], TimeNames, []);
                FilterImage = [];
                TimePoints = [];
            end
            if ~NOPLOT
                gcfpos = [120    50   760   700];
                set(gcf,'Position',gcfpos);
                drawnow
            end
        else
            if ~isempty(FileName)
                [FilterImage, TimePoints] = generate_filter_image(Fmat,Classname,FmatNchan,[],thres,FmatChanNames,TimeNames,FileName,ClassifierDesc,SWsortchans,CommentTextCell,NOPLOT,[]);
            else
                [FilterImage, TimePoints] = generate_filter_image(Fmat,Classname,FmatNchan,[],thres,FmatChanNames,TimeNames,[],ClassifierDesc,SWsortchans,CommentTextCell,NOPLOT,[]);
            end
            if ~NOPLOT
                gcfpos = [100    50   760   700];
                set(gcf,'Position',gcfpos);
                drawnow
            end
        end
    end
end



return