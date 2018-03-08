function yt = getyticks (YVal, ystart, yint, yend)
if nargin == 2
    ylist = ystart;
else
    ylist = ystart:yint:yend;
end
yt = zeros(1,length(ylist));
for i = 1:length(ylist)
    [~,yt(i)] = min(abs(ylist(i) - YVal));
end
yt = unique(yt);

