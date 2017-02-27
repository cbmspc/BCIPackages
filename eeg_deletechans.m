function [data, ChanNames] = eeg_deletechans (data, ChanNames, DelChanNames)
% Deletes channel(s) from data(chan,time) and ChanNames(cell)
ChanNames = upper(ChanNames);
DelChanNames = upper(DelChanNames);

if isstr(ChanNames)
    ChanNames = string_to_cell(ChanNames,' ,-_');
end
if isstr(DelChanNames)
    DelChanNames = string_to_cell(DelChanNames,' ,-_');
end

for i = 1:length(DelChanNames)
    EName = DelChanNames{i};
    Nchan = length(ChanNames);
    idxall = 1:Nchan;
    idxdel = find(strcmpi(ChanNames,EName));
    idxkeep = setdiff(idxall,idxdel);
    ChanNames = ChanNames(idxkeep);
    data = data(idxkeep,:);
end

