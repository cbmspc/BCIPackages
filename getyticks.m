function yt = getyticks (YVal, ystart, yint, yend, Interpolate)
if nargin == 2 || length(ystart) > 1
    ylist = ystart;
else
    ylist = ystart:yint:yend;
end

if ~exist('Interpolate','var') || isempty(Interpolate)
    Interpolate = 0;
end

if Interpolate
    yt = interp1(YVal, 1:length(YVal), ylist, 'linear', 'extrap');
else
    yt = zeros(1,length(ylist));
    for i = 1:length(ylist)
        [~,yt(i)] = min(abs(ylist(i) - YVal));
    end
end
yt = unique(yt);
