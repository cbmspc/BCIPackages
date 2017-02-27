function idx = eeg_findcoord (CoordList, Coords)
% Find the index numbers of (x,y) coordinates in the CoordList list
NC = size(Coords,1);
for i = 1:NC
    idx_xmatch = find(CoordList(:,1) == Coords(i,1));
    idx_ymatch = find(CoordList(:,2) == Coords(i,2));
    o = intersect(idx_xmatch,idx_ymatch);
    if length(o) ~= 1
        o = NaN;
    end
    idx(i) = o;
end
