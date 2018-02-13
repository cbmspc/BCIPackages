function yt = getyticks (YVal, ystart, yint, yend)
ylist = ystart:yint:yend;
yt = zeros(1,length(ylist));
for i = 1:length(ylist)
    [~,yt(i)] = min(abs(ylist(i) - YVal));
end
yt = unique(yt);

