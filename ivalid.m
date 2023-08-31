function [pcorrect, pconfuse] = ivalid (TrainData, TrainLabels, TestData, TestLabels, varargin)

OPTS = 1;
DRFUN = 'cpca';
DRPARM = {[]};
FEFUN = {'aida','aida','aida','aida','aida','lda'};
FEPARM = {5,4,3,2,1,1};
CFUN = {'',''};
CNAME = {'quadratic','linear'};
PRIOR = 'equal';
CPARM = {[],[]};
SWsilent = 0;
sdfactor = 1;
COLUMN = 2;
FPID = 1;
SWlogscale = 0;
BALANCED = 1;
if ~iscell(DRFUN)
    DRFUN = {DRFUN};
end

TrainLabels = TrainLabels(:);
TestLabels = TestLabels(:);

if mod(length(varargin),2)
    error('Properties list must be in pairs, ie. property1 name, property1 value, ...');
end

for i = 1:2:length(varargin)
    switch upper(varargin{i})
        case 'DRFUN'
            DRFUN = lower(varargin{i+1});
        case 'DRPARM'
            if iscell(varargin{i+1}) && ~isempty(varargin{i+1}) && iscell(varargin{i+1}{1})
                % 20150316 Double-wrap fix
                DRPARM = varargin{i+1}{1};
            else
                DRPARM = varargin{i+1};
            end
        case 'PRIOR'
            PRIOR = varargin{i+1};
        case 'CPARM'
            CPARM = varargin{i+1};
        case 'TIME'
            TIMESPEC = varargin{i+1};
        case 'TIMESPEC'
            TIMESPEC = varargin{i+1};
        case 'CHAN'
            CHANSPEC = varargin{i+1};
        case 'CHANSPEC'
            CHANSPEC = varargin{i+1};
        case 'LABELSPEC'
            LABELSPEC = varargin{i+1};
        case 'OPTS'
            OPTS = varargin{i+1};
        case 'FPRINTID'
            FPID = varargin{i+1};
        case 'LOGSCALE'
            if varargin{i+1} == 1
                SWlogscale = 1;
            elseif varargin{i+1} == 0
                SWlogscale = 0;
            end
        case 'BALANCED'
            if varargin{i+1} == 2
                BALANCED = 2;
            elseif varargin{i+1} == 1
                BALANCED = 1;
            elseif varargin{i+1} == 0
                BALANCED = 0;
            end
        case 'SILENT'
            if varargin{i+1} == 1
                SWsilent = 1;
                OPTS = 0;
            end
        case 'COLUMN'
            if varargin{i+1} == 2
                COLUMN = 2;
            elseif varargin{i+1} == 1
                COLUMN = 1;
            end
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


if (isempty(who('TestLabels')) || isempty(TestLabels)) && ~isempty(TestData) && isstruct(TestData(1))
    if OPTS > 0
        fprintf(FPID, 'Reconstructing rawdata... ');
    end
    [TestData,TestLabels,ChanNames,samplerate,datevector] = bsmspack_convert_to_rawdata(TestData,'center');
    [S T N] = size(TestData);
    if OPTS > 0
        fprintf(FPID, ' %i sensors, %i time points, %i trials.\n',S,T,N);
    end
else
    if size(TestData,3) > 1
        [S T N] = size(TestData);
    else
        [N T] = size(TestData);
        S = 1;
    end
end


Data = [];
Labels = [];

if size(TrainData,3) == 1
    TrainDataDim = size(TrainData,2);
else
    TrainDataDim = size(TrainData,1) * size(TrainData,2);
end

if size(TestData,3) == 1
    TestDataDim = size(TestData,2);
else
    TestDataDim = size(TestData,1) * size(TestData,2);
end

if TrainDataDim ~= TestDataDim
    error('Dimension mismatch between TrainData and TestData');
end


if ~isempty(who('TIMESPEC')) || ~isempty(who('CHANSPEC'))
   if size(TrainData,3) == 1 || size(TestData,3) == 1
       warning('TrainData or TestData is already reshaped. TIMESPEC and CHANSPEC are ignored!');
   else
       if isempty(who('TIMESPEC'))
           TIMESPEC = [1:size(TrainData,2)];
       end
       if isempty(who('CHANSPEC'))
           CHANSPEC = [1:size(TrainData,1)];
       end
       TrainData = TrainData(CHANSPEC,TIMESPEC,:);
       TestData = TestData(CHANSPEC,TIMESPEC,:);
       if OPTS > 0
           disp('Subsets of TrainData and TestData extracted using:');
           disp(['CHANSPEC = ' summarize_idx(CHANSPEC)]);
           disp(['TIMESPEC = ' summarize_idx(TIMESPEC)]);
       end
   end
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
    
    if any(TestData(:) <= 0)
        TestData = eps+abs(TestData);
        warning('TestData contains non-positive numbers. To compute log, it has been fixed with abs and eps.');
    end
    TestData = log(TestData);
    if OPTS > 0
        disp('TestData time/freq series has been converted to natural-log scale.');
    end
end



% Automatically delete data that do not have corresponding labels in the
% other set of data (i.e., delete the odd labels)
oddlabels = setxor(unique(TestLabels),unique(TrainLabels));
if ~isempty(oddlabels)
    for i = 1:length(oddlabels)
        tsidxx{i} = find(TestLabels==oddlabels(i));
        tridxx{i} = find(TrainLabels==oddlabels(i));
    end
    atsidxx = cat(1,tsidxx{:});
    atridxx = cat(1,tridxx{:});
    TestLabels = TestLabels(setdiff(1:length(TestLabels),atsidxx));
    TrainLabels = TrainLabels(setdiff(1:length(TrainLabels),atridxx));
    if size(TestData,3) == 1
        TestData = TestData(setdiff(1:size(TestData,1),atsidxx),:);
    else
        TestData = TestData(:,:,setdiff(1:size(TestData,3),atsidxx));
    end
    if size(TrainData,3) == 1
        TrainData = TrainData(setdiff(1:size(TrainData,1),atridxx),:);
    else
        TrainData = TrainData(:,:,setdiff(1:size(TrainData,3),atridxx));
    end
    if OPTS > 0
        disp('Deleted test and train data with non-mutual labels.');
        for i = 1:length(oddlabels)
            disp(['LABEL = ' num2str(oddlabels(i)) ' (DELETED)']);
        end
    end
    clear atsidxx tsidxx atridxx tridxx
end

% Make sure LABELSPEC only contains labels in TrainLabels. As a result of
% the previous block of code, this will also imply making sure LABELSPEC
% only contains labels in TestLabels as well.
if ~isempty(who('LABELSPEC'))
    LABELSPEC = unique(LABELSPEC(:));
    oddlabelspec = setdiff(LABELSPEC,unique(TrainLabels));
    if ~isempty(oddlabelspec)
        for i = 1:length(oddlabelspec)
            lsidxx{i} = find(LABELSPEC==oddlabelspec(i));
        end
        alsidxx = cat(1,lsidxx{:});
        LABELSPEC = LABELSPEC(setdiff(1:length(LABELSPEC),alsidxx));
        if OPTS > 0
            disp('Deleted labels in LABELSPEC with no corresponding labels in train data.');
            for i = 1:length(oddlabelspec)
                disp(['LABELSPEC = ' num2str(oddlabelspec(i)) ' (DELETED)']);
            end
        end
        clear alsidxx lsidxx
        
    end
end

if ~isempty(who('LABELSPEC'))
    for i = 1:length(LABELSPEC)
        tridx{i} = find(TrainLabels==LABELSPEC(i));
        tsidx{i} = find(TestLabels==LABELSPEC(i));
    end
    atridx = cat(1,tridx{:});
    atsidx = cat(1,tsidx{:});
    TrainLabels = TrainLabels(atridx);
    TestLabels = TestLabels(atsidx);
    if size(TrainData,3) == 1
        TrainData = TrainData(atridx,:);
    else
        TrainData = TrainData(:,:,atridx);
    end
    if size(TestData,3) == 1
        TestData = TestData(atsidx,:);
    else
        TestData = TestData(:,:,atsidx);
    end
    if OPTS > 0
        disp('Only the following classes are included:');
        for i = 1:length(LABELSPEC)
            disp(['LABEL = ' num2str(LABELSPEC(i)) ' (' num2str(length(tridx{i})) ' training observations, ' num2str(length(tsidx{i})) ' test observations)']);
        end
    end
end

if size(TrainData,3) > 1
    if OPTS > 0
        disp('Reshaping Training data...');
    end
    TrainData = rdreshape(TrainData);
end

if size(TestData,3) > 1
    if OPTS > 0
        disp('Reshaping Test data...');
    end
    TestData = rdreshape(TestData);
end


% Balancing
clear pcorrectC pconfuseC pcorrect pconfuse
if BALANCED
    MRUN = 10;
    if OPTS > 0
        disp('Randomly subsampling the TrainData so that each class has the same number of trials.');
        disp('TestData is not subsampled. All TestData are included.');
    end
    classes = unique(TrainLabels);
    Nclass = length(classes);
    for j = Nclass:-1:1
        NtrialTrain(j) = length(find(TrainLabels==classes(j)));
    end
    NtrialTrainMin = min(NtrialTrain);
    
    %wh = waitbar(0, 'Inter validating');
    for i = 1:MRUN
        clear TrainDataC TrainLabelsC
        for j = Nclass:-1:1
            % First convert to cell
            TrainDataC{j} = TrainData(TrainLabels==classes(j),:);
            TrainLabelsC{j} = classes(j)*ones(NtrialTrain(j),1);
            if NtrialTrain(j) > NtrialTrainMin
                % Subsample this class.
                keepidx = randperm(NtrialTrain(j),NtrialTrainMin);
                TrainDataC{j} = TrainDataC{j}(keepidx,:);
                TrainLabelsC{j} = TrainLabelsC{j}(keepidx,:);
            end
        end
        % Then unconvert from cell
        TrainDataSub = cat(1,TrainDataC{:});
        % All labels are also rearranged
        TrainLabelsSub = cat(1,TrainLabelsC{:});
        
        Data = cat(1,TestData,TrainDataSub);
        Labels = cat(1,TestLabels,TrainLabelsSub);
        TestIdx = [1:length(TestLabels)];
        
        %Po 2023-07-10: If command fails (e.g. Inf/NaN), skip it for the
        %averaging
        try
            [pcorrectC{i},pconfuseC{i}] = dataproc_main_multiintervalidation(Data,Labels,TestIdx,DRFUN,DRPARM,FEFUN,FEPARM,CFUN,CNAME,PRIOR,CPARM,OPTS);
        catch
            pcorrectC{i} = NaN;
            pconfuseC{i} = NaN;
        end
        
        %if ishandle(wh)
        %    try
        %        waitbar(i/MRUN, wh, sprintf('Inter validating, last Pcorrect=%-5.3f', median(median(cell2mat(pcorrectC{i})))));
        %    end
        %end
    end
    % 2021-10-20 Lazy: Combine them and forget it
    
    RealMRUN = MRUN;
    for i = 1:MRUN
        if ~iscell(pcorrectC{i}) && isnan(pcorrectC{i})
            RealMRUN = RealMRUN - 1;
        end
    end
    jmax = numel(DRPARM);
    kmax = numel(FEPARM);
    lmax = numel(CPARM);
    if RealMRUN > 0
        for i = 1:MRUN
            if iscell(pcorrectC{i})
                for j = size(pcorrectC{i},1):-1:1
                    for k = size(pcorrectC{i},2):-1:1
                        for l = size(pcorrectC{i},3):-1:1
                            if i == 1
                                pcorrect{j,k,l} = pcorrectC{i}{j,k,l} / RealMRUN;
                                pconfuse{j,k,l} = pconfuseC{i}{j,k,l} ./ RealMRUN;
                            else
                                pcorrect{j,k,l} = pcorrect{j,k,l} + pcorrectC{i}{j,k,l} / RealMRUN;
                                pconfuse{j,k,l} = pconfuse{j,k,l} + pconfuseC{i}{j,k,l} / RealMRUN;
                            end
                        end
                    end
                end
            end
        end
    else
        pcorrect = cell(jmax,kmax,lmax);
        pcorrect{:} = 0;
        pconfuse = cell(jmax,kmax,lmax);
        pcorrect{:} = (double(~eye(Nclass)))/(Nclass-1);
    end
    
    %if ishandle(wh)
    %    try
    %        delete(wh);
    %        drawnow
    %    end
    %end
    
    
else
    TrainDataSub = TrainData;
    TrainLabelsSub = TrainLabels;
    
    Data = cat(1,TestData,TrainDataSub);
    Labels = cat(1,TestLabels,TrainLabelsSub);
    TestIdx = [1:length(TestLabels)];
    
    
    [pcorrect,pconfuse] = dataproc_main_multiintervalidation(Data,Labels,TestIdx,DRFUN,DRPARM,FEFUN,FEPARM,CFUN,CNAME,PRIOR,CPARM,OPTS);

    
end









classes = unique(Labels);
Nclass = length(classes);

for j = 1:Nclass
    NtrialTest(j) = length(find(TestLabels==classes(j)));
    NtrialTrain(j) = length(find(TrainLabels==classes(j)));
end

if ~SWsilent
    fprintf(FPID, '\n');
    if strcmpi(DRFUN,'cpca')
        cpcamat1 = dataproc_func_cpca(TrainData,TrainLabels);
        for i = 1:length(cpcamat1)
            cpcadim1(i) = size(cpcamat1{i},2);
        end
        cpcamat2 = dataproc_func_cpca(TestData,TestLabels);
        for i = 1:length(cpcamat2)
            cpcadim2(i) = size(cpcamat2{i},2);
        end
    elseif strcmpi(DRFUN,'pca')
        pcamat1 = dataproc_func_pca(TrainData,TrainLabels);
        pcadim1 = size(pcamat1,2);
        pcamat2 = dataproc_func_pca(TestData,TestLabels);
        pcadim2 = size(pcamat2,2);
    end
    fprintf(FPID, 'Inter Validation Report\n');
    
    fprintf(FPID, '    Input dimension: %i\n',size(TrainData,2));
    %fprintf(FPID, 'Dimension reduction: %s (%g)\n',DRFUN{1},DRPARM{1});
    if iscell(DRFUN)
        fprintf(FPID, 'Dimension reduction: %s (%g)\n',DRFUN{1},DRPARM{1});
    else
        fprintf(FPID, 'Dimension reduction: %s (%g)\n',DRFUN,DRPARM);
    end

    fprintf(FPID, 'TRAINING\n');
    fprintf(FPID, '   Number of trials: ');
    for i = 1:Nclass-1
        fprintf(FPID, '%i (class %i), ',NtrialTrain(i),classes(i));
    end
    fprintf(FPID, '%i (class %i)\n',NtrialTrain(Nclass),classes(Nclass));
    if strcmpi(DRFUN,'cpca')
        fprintf(FPID, ' Reduced dimensions: ');
        for i = 1:length(cpcadim1)-1
            fprintf(FPID, '%i (class %i), ',cpcadim1(i),classes(i));
        end
        fprintf(FPID, '%i (class %i)\n',cpcadim1(end),classes(end));
    elseif strcmpi(DRFUN,'pca')
        fprintf(FPID, '  Reduced dimension: ');
        fprintf(FPID, '%i\n',pcadim1);
    end
    fprintf(FPID, 'TEST\n');
    fprintf(FPID, '   Number of trials: ');
    for i = 1:Nclass-1
        fprintf(FPID, '%i (class %i), ',NtrialTest(i),classes(i));
    end
    fprintf(FPID, '%i (class %i)\n',NtrialTest(Nclass),classes(Nclass));
    if strcmpi(DRFUN,'cpca')
        fprintf(FPID, ' Reduced dimensions: ');
        for i = 1:length(cpcadim2)-1
            fprintf(FPID, '%i (class %i), ',cpcadim2(i),classes(i));
        end
        fprintf(FPID, '%i (class %i)\n',cpcadim2(end),classes(end));
    elseif strcmpi(DRFUN,'pca')
        fprintf(FPID, '  Reduced dimension: ');
        fprintf(FPID, '%i\n',pcadim2);
    end
    fprintf(FPID, 'Prior probabilities: [%s]\n',num2str(PRIOR));
end

pcor = shiftdim(pcorrect(1,:,:),1);
pcon = shiftdim(pconfuse(1,:,:),1);

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

%[Y,I] = sort(max(pcormean,[],2));
%clear('Y');
%I = flipud(I(:));

[Y,I] = sort(reshape(pcormean,[],1),1,'descend');
[Y J] = sort(I,1,'ascend');
J = reshape(J,[],2);

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
                    fprintf(FPID, '%5.1f%%        %5s',pcormean(fi,cf)*100,RANKTEXT);
                    %disp([' FE: ' FEFUN{fi} '(' num2str(FEPARM{fi}) ')' ' ' 'CF: ' CNAME{cf} ': Pcorrect = ' num2str(pcormean(fi,cf)*100,'%5.1f') '%' ' ' RANKTEXT]);
                for i = 3:Nclass
                    fprintf(FPID, '               ');
                end
                fprintf(FPID, ' | ');
                %disp('Confusion matrix');
            end
            fprintf(FPID, '\n');
            me = pconmean(:,:,fi,1)*100;
            sd = pconstd(:,:,fi,1)*100;
            for i = 1:size(me,1)
                for cf = 1:length(CNAME)
                    me = pconmean(:,:,fi,cf)*100;
                    sd = pconstd(:,:,fi,cf)*100;
                    for j = 1:size(me,2)
                            fprintf(FPID, '%5.1f%%         ', me(i,j));
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
                    fprintf(FPID, 'Pcorrect =%5.1f%% %s\n',pcormean(fi,cf)*100,RANKTEXT);
                    %disp([' FE: ' FEFUN{fi} '(' num2str(FEPARM{fi}) ')' ' ' 'CF: ' CNAME{cf} ': Pcorrect = ' num2str(pcormean(fi,cf)*100,'%5.1f') '%' ' ' RANKTEXT]);
                %disp('Confusion matrix');
                me = pconmean(:,:,fi,cf)*100;
                sd = pconstd(:,:,fi,cf)*100;
                for i = 1:size(me,1)
                    for j = 1:size(me,2)
                            fprintf(FPID, '%5.1f%%', me(i,j));
                    end
                    fprintf(FPID, '\n');
                end
                fprintf(FPID, '\n');
            end
        end
    end
    fprintf(FPID, '\n');
end
