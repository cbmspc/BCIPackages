% Find the envelope indexes of individual pulse trains (3rd column of SmallerBounds)
% Also find the number of pulses in each envelope (3rd column of BiggerBounds)
% Bounds are boundaries, i.e. [start end] pairs each line
%
% Example: 
%  bbs = getbounds2(pa(:,1), Fs, 'ManualThreshold', 2, 'MergeLen', 1.00, 'SWloadManualCoords', 1, 'Squared', -1, 'Fcutoff', [0 0], 'MinLen', 0.001, 'SWskipmanualinstruction', 1, 'SWautoaccept', 0);
%  bbs_tiny = getbounds2(pa(:,1), Fs, 'ManualThreshold', 2, 'MergeLen', 0.1, 'SWloadManualCoords', 1, 'Squared', -1, 'Fcutoff', [0 0], 'MinLen', 0.001, 'SWskipmanualinstruction', 1, 'SWautoaccept', 0);
%  [bbs_tiny, bbs] = findboundsinbounds(bbs_tiny, bbs);
%  bbs = bbs(bbs(:,3) >= 5,:); % only keep envelopes with at least 5 pulses
%


function [SmallerBounds, BiggerBounds] = findboundsinbounds (SmallerBounds, BiggerBounds)

BiggerBounds = uniquerows(BiggerBounds);

BiggerBounds(:,3) = 0;

for i = 1:size(SmallerBounds,1)
    l = SmallerBounds(i,1) > BiggerBounds(:,1) & SmallerBounds(i,2) < BiggerBounds(:,2);
    if nnz(l)
        ind = find(l,1);
        SmallerBounds(i,3) = ind;
        BiggerBounds(ind,3) = BiggerBounds(ind,3) + 1;
    else
        SmallerBounds(i,3) = 0;
    end
    
end