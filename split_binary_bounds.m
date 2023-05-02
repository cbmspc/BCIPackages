% Split experimental data into multiple trials for each boundary
% For example, if you recorded an experiment with alternating 10 seconds
% idle and move epochs, you can split each 10-sec epoch into:
% 1. Ignore the first 1 second to account for delay
% 2. Split the remaining 9 seconds into 3-second trials.
% The durations of epochs do not have to be the same.
%
function [BoundsSec, Labels] = split_binary_bounds (IdleBoundsSec, MoveBoundsSec, SkipFirstSec, NSplitsWithinABound)
Ni = size(IdleBoundsSec,1);
Nm = size(MoveBoundsSec,1);
BoundsSec = nan(NSplitsWithinABound*(Ni+Nm),2);
Labels = nan(size(BoundsSec,1),1);
k = 0;
for i = 1:Ni
    a = IdleBoundsSec(i,1) + SkipFirstSec;
    b = IdleBoundsSec(i,2);
    dur = (b-a)/NSplitsWithinABound;
    for j = 0:NSplitsWithinABound-1
        k = k + 1;
        BoundsSec(k,:) = [a+dur*j, a+dur*(j+1)];
        Labels(k) = 0;
    end
end
for i = 1:Nm
    a = MoveBoundsSec(i,1) + SkipFirstSec;
    b = MoveBoundsSec(i,2);
    dur = (b-a)/NSplitsWithinABound;
    for j = 0:NSplitsWithinABound-1
        k = k + 1;
        BoundsSec(k,:) = [a+dur*j, a+dur*(j+1)];
        Labels(k) = 1;
    end
end
