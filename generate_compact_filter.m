function SFilter = generate_compact_filter (trfmat, EEGChanNames, ChanIdx, TimeTicks, FileName)
SWstaggerchannames = 1;
Nchan = length(EEGChanNames);
if Nchan > 34
    SWstaggerchannames = 1;
end
Nclass = length(trfmat);
if ~exist('ChanIdx','var') || isempty(ChanIdx)
    ChanIdx = 1:Nchan;
end
for i = 1:Nclass
    Filter{i} = reshape(trfmat{i},Nchan,[]);
end
Ntime = size(Filter{1},2);

for C = 1:2^(Nclass-1)
    V = getbinarycombnum(Nclass,C);
    SFilter = zeros(Nchan,Ntime);
    for i = 1:Nclass
        SFilter = SFilter + V(i)*Filter{i};
    end
    SFN(C) = norm(SFilter);
end

[Y C] = max(SFN);
V = getbinarycombnum(Nclass,C);
SFilter = zeros(Nchan,Ntime);
for i = 1:Nclass
    SFilter = SFilter + V(i)*Filter{i};
end

SFilter = SFilter ./ max(abs(SFilter(:)));

MaxChanNameLength = -1;
for i = 1:length(EEGChanNames)
    MaxChanNameLength = max(MaxChanNameLength,length(EEGChanNames{i}));
end

DisplayedEEGChanNames = EEGChanNames;
if SWstaggerchannames
    for i = 2:2:length(EEGChanNames)
        
        if i < Nchan
        nn = max(length(DisplayedEEGChanNames{i-1}),length(DisplayedEEGChanNames{i+1}));
        else
            nn = length(DisplayedEEGChanNames{i-1});
        end
        
        DisplayedEEGChanNames{i} = [DisplayedEEGChanNames{i} ' ' repmat(' ',1,nn)];
    end
end

if nargout == 0
    figure(10029);
    imagesc(SFilter(ChanIdx,:),[-1 1]);
    NIchan = length(ChanIdx);
    colorbar;
    set(gca,'YTick',[1:NIchan],'YTickLabel',DisplayedEEGChanNames(ChanIdx));
    if SWstaggerchannames
        set(gca, 'FontName', 'Consolas', 'FontSize', max(12,min(18,round(18*64/Nchan))));
    else
        set(gca, 'FontName', 'Consolas', 'FontSize', max(6,min(16,round(12*64/Nchan))));
    end
    
    %set(gca,'XTick',[0:Ntime/8:Ntime]);
    %dT = (Time2-Time1)/8;
    %set(gca,'XTickLabel',[Time1:dT:Time2])
    set(gca,'XTick',[1:length(TimeTicks)]);
    set(gca,'XTickLabel',TimeTicks);
end

if ~isempty(who('FileName')) && ~isempty(FileName)
    saveas(gcf, [FileName ' compact filter image'], 'fig');
    if ispc
        print(gcf,'-dpng',[FileName ' compact filter image']);
    end
end







function V = getbinarycombnum (N, C)
S = dec2bin(C);
for i = 1:(N-length(S))
    S = ['0' S];
end
for i = 1:length(S)
    switch S(i)
        case '0'
            V(i) = -1;
        case '1'
            V(i) = 1;
    end
end
