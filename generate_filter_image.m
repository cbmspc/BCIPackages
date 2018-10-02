%

function [FilterImage, XTickLabel] = generate_filter_image (TrainData, TrainLabels, Nchan, Nm, threshold, ChanNames, TimeNames, FileName, ClassifierDesc, SWsortchans, CommentTextCell, SWnoplot, FigPos)
% Nm = number of final dimensions
% threshold = number of standard deviations
% ChanNames = cell array of channel names

if isempty(who('Nm')) || isempty(Nm)
    Nm = 1;
end

if isempty(who('threshold')) || isempty(threshold)
    threshold = -inf;
end

if ~exist('SWsortchans','var') || isempty(SWsortchans)
    SWsortchans = 0;
end

if length(ChanNames) ~= Nchan
    SWsortchans = 0;
end

if ~exist('FigPos', 'var') || length(FigPos) ~= 4
    FigPos = [];
end

if ~isempty(who('ChanNames')) && ~isempty(ChanNames)
    SWlabelchannels = 1;
    % Ensure it is horizontal
    ChanNames = ChanNames(:).';
    
    % Sort channels according to a master list
    if SWsortchans
        ChanNames_sorted = getsortedchannames();
        ChanNames_sorted = ChanNames_sorted(:).';
        
        % Empty channels are not shown
        tmp = chan2idx(ChanNames, ChanNames_sorted);
        MasterChanNames = ChanNames(tmp(tmp>0));
        
        
        % Add channels that only exists in ChanNames (input argument) to end
        MasterChanNames = [MasterChanNames ChanNames(chan2idx(MasterChanNames, ChanNames)==0)];
        
        % The channels on the final image will be according to this
        MasterChanIdx = chan2idx(ChanNames, MasterChanNames);
        
    end
else
    SWlabelchannels = 0;
end

if ~exist('TimeNames', 'var')
    TimeNames = {};
end

if isempty(TrainLabels) || iscell(TrainLabels) % TrainLabels doubles as classname if it is a cell
    FilterImage = TrainData;
    Ns = size(FilterImage,2);
    Nchan = size(FilterImage{1},1);
else

    if size(TrainData,3) > 1
        disp('Reshaping data...');
        Nchan = size(TrainData,1);
        TrainData = rdreshape(TrainData);
    end


    trfmat = cpca_aida(TrainData,TrainLabels,Nm);
    
    ClassifierDesc = 'CPCA+AIDA';
    
%     cpcamat = dataproc_func_cpca(TrainData,TrainLabels);
%     ldamat{1} = dataproc_func_lda(TrainData*cpcamat{1},TrainLabels,1);
%     ldamat{2} = dataproc_func_lda(TrainData*cpcamat{2},TrainLabels,1);
%     trfmat{1} = cpcamat{1}*ldamat{1};
%     trfmat{2} = cpcamat{2}*ldamat{2};
    
    Ns = 1;
    if iscell(trfmat)
        Ns = length(trfmat);
    end
end

if ~exist('ClassifierDesc','var')
    ClassifierDesc = '';
end

Negated = zeros(Nm,Ns);

for m = 1:Nm
    for s = 1:Ns
        if ~isempty(TrainLabels) && ~iscell(TrainLabels)
            FilterImage{m,s} = reshape(trfmat{s}(:,m),Nchan,[]);
        end
        AF{m,s} = abs(FilterImage{m,s});
        RAF{m,s} = reshape(AF{m,s},1,[]);
        mm(m,s) = mean(RAF{m,s});
        ss(m,s) = std(RAF{m,s});
        thres(m,s) = mm(m,s) + threshold * ss(m,s);
        BF{m,s} = FilterImage{m,s};
        BF{m,s}(find(abs(FilterImage{m,s}) < thres(m,s))) = 0;
        cmin(m,s) = min(min(BF{m,s}));
        cmax(m,s) = max(max(BF{m,s}));
    end
    for s = 2:Ns
        if norm(BF{m,1} + BF{m,s}) < norm(BF{m,1} - BF{m,s})
            % If image difference is smaller for negated image
            BF{m,s} = -BF{m,s};
            Negated(m,s) = 1;
        end
        cmin(m,s) = min(min(BF{m,s}));
        cmax(m,s) = max(max(BF{m,s}));
    end

end

cmin = min(cmin,[],2);
cmax = max(cmax,[],2);

% Normalize displayed pixels to [-1,1]
NBF = BF;
for m = 1:Nm
    for s = 1:Ns
        NBF{m,s} = NBF{m,s} ./ max(abs(NBF{m,s}(:)));
    end
end

close all;

%AxisWiggleRoomLR = [0.065 0.01];
%AxisWiggleRoomDU = [0.03 0.03];
%AxisWiggleRoomLR = [0.075 0.03];
AxisWiggleRoomLR = [0.105 0.03];
AxisWiggleRoomDU = [0.03 0.03];


for s = Ns:-1:1
    AxisPos(s,:) = [AxisWiggleRoomLR(1), (Ns-s)*(1.0/Ns)+AxisWiggleRoomDU(1), 1.0-sum(AxisWiggleRoomLR), 1.0/Ns-sum(AxisWiggleRoomDU)];
end

if SWlabelchannels
    if SWsortchans
        % Re-arrange data matrix according to master channel list
        New_Nchan = length(MasterChanNames);
        for m = 1:Nm
            for s = 1:Ns
                tmp = zeros(New_Nchan, size(NBF{m,s},2));
                for ch = 1:New_Nchan
                    if MasterChanIdx(ch) > 0
                        tmp(ch,:) = NBF{m,s}(MasterChanIdx(ch),:);
                    end
                end
                NBF{m,s} = tmp;
            end
        end
        Nchan = New_Nchan;
        ChanNames = MasterChanNames;
    end
    
    % Stagger channel names for easier reading
    ChanNameLen = [];
    for i = 1:length(ChanNames)
        ChanNameLen(i) = length(ChanNames{i});
    end
    %stagger = 1 + ceil(max(ChanNameLen)*2.2);
    stagfac = 2.1;
    
    stagger = 1 + ceil(max(ChanNameLen)*stagfac);
    for i = 1:2:length(ChanNames)
        ChanNames{i} = [ChanNames{i} char(ones(1,stagger)*' ')];
    end
end


if ~exist('SWnoplot','var') || isempty(SWnoplot)
    SWnoplot = 0;
end

% Actual plotting
fighand = zeros(1,Nm);
for m = 1:Nm
    fighand(m) = figure(m);
    clf
    set(fighand(m),'Position',[100+(m-1)*100 50 800 768]);
    if SWnoplot
        set(fighand(m), 'Visible', 'off');
    end
    for s = 1:Ns
        %subplot(Ns,1,s);
        shand(m,s) = subplot('Position',AxisPos(s,:));
        imagesc(real(NBF{m,s}),[-1,1]);
        %imagesc(BF{m,s});
        %caxis([cmin(m) cmax(m)]);
        
        XTickLabel = {};
        if ~isempty(TimeNames)
            if numel(TimeNames) == size(BF{m,s},2)
                % Exact points are specified
                XTicks = 1:size(BF{m,s},2);
                XTickLabel = TimeNames;
                %20111018 Changed 21 to 1001
                if length(XTicks) > 1001
                    keeptick = round(linspace(1,length(XTicks),21));
                    XTicks = XTicks(keeptick);
                    XTickLabel = XTickLabel(keeptick);
                end
                set(gca,'XTick',XTicks,'XTickLabel',XTickLabel);
            elseif isnumeric(TimeNames) && numel(TimeNames) == 2
                % Left and right bounds are specified
                Ntp = size(BF{m,s},2);
                %TimeOffset = TimeNames(1);
                %TimeFactor = (TimeNames(2)-TimeOffset)/Ntp;
                TimeDiv = diff(TimeNames)/10;
                XTicks = linspace(0,Ntp,11);
                %XTicks(1) = 1;
                XTickLabel = TimeNames(1):TimeDiv:TimeNames(2);
                set(gca,'XTick',XTicks,'XTickLabel',XTickLabel);
                %set(gca,'XTick',XTicks,'XTickLabel',XTicks*TimeFactor+TimeOffset);
            else
                % Neither
                warning('Dimension mismatch between TimeNames and Filter');
            end
        end
        
        %xlabel('Time or Sample Points');
        Sname = num2str(s-1);
        if ~isempty(TrainLabels)
            if ~iscell(TrainLabels)
                classes = unique(TrainLabels);
                Sname = ['"' num2str(classes(s)) '"'];
            else
                Sname = ['"' num2str(TrainLabels{s}) '"'];
            end
        end
        %ylabel(['Class Subspace ' Sname], 'FontName', 'Consolas', 'FontSize', 10);
        ylabel(['Class Subspace ' Sname]);

        if SWlabelchannels
            set(gca,'YTick',(1:1:Nchan),'YTickLabel',ChanNames(1:1:Nchan));
            %set(gca,'FontSize',min([16,floor(48/sqrt(Nchan))]));
        end
        set(gca,'FontSize',max(4,min([10,floor(48/max(sqrt(Nchan),2*sqrt(size(BF{m,s},2))))])));
        %set(gca, 'FontSize', 10, 'FontName', 'Consolas');

        optext1 = '';
        if Negated(m,s)
            optext1 = ' (negated)';
        end

        optext2 = '';
        if threshold > -inf
            optext2 = [', thres=\mu+' num2str(threshold) '\sigma'];
        else
            optext2 = [', no thres'];
        end
        
        optext3 = '';
        if exist('CommentTextCell','var') && ~isempty(CommentTextCell)
            if ~iscell(CommentTextCell)
                tmp = cell(1,Ns);
                tmp(1:Ns) = {CommentTextCell};
                CommentTextCell = tmp;
            end
            optext3 = [', ' CommentTextCell{s}];
        end

        %title([ClassifierDesc ' Filter' optext1 ', subspace ' Sname ', dim ' num2str(m) ' of ' num2str(Nm) optext2 optext3], 'FontSize', 9);
        %title([ClassifierDesc ' Filter' optext1 ', dim ' num2str(m) ' of ' num2str(Nm) optext2 optext3], 'FontSize', 10, 'FontName', 'Consolas');
        title([ClassifierDesc ' Filter' optext1 ', dim ' num2str(m) ' of ' num2str(Nm) optext2 optext3], 'FontSize', 9);
        cb = colorbar('EastOutside','FontSize',8);
        set(cb, 'YTick', [-1:0.1:1]);
        %if SWlabelchannels
        %   ax2 = axes('Position', get(gca, 'Position'), 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'Color', 'none', 'XColor', 'k', 'YColor', 'k', 'XTick', [], 'YDir', 'reverse', 'YLim', get(gca,'YLim'), 'FontSize', min([16,floor(48/sqrt(Nchan))]));
        %   set(ax2,'YTick',[2:2:Nchan],'YTickLabel',ChanNames(2:2:Nchan));
        %end
    end
end




if length(FilterImage) == 1
    G = FilterImage{1};
    clear('FilterImage');
    FilterImage = G;
end

if ~isempty(FigPos)
    for m = 1:Nm
        if ishandle(fighand(m))
            set(fighand(m), 'Position', FigPos);
            pause(0.1);
        end
    end
end

if ~isempty(who('FileName')) && ~isempty(FileName)
    Ne = length(dir([FileName ' (m=*)*.png']));
    Fsuffix = '';
    if Ne > 0
        Fsuffix = [' [' num2str(Ne+1) ']'];
    end
    for m = 1:Nm
        saveas(m, [FileName ' (m=' num2str(m) ')' Fsuffix '.fig'], 'fig');
        %if ispc
        %   print(m,'-dpng',[FileName ' (m=' num2str(m) ')' Fsuffix '.png']);
        %end
    end
end

if SWnoplot
    for m = 1:Nm
        if ishandle(fighand(m))
            delete(fighand(m));
        end
    end
end

