function [yt, intervals] = colorbarticks (CLim, Div)
% The default colorbar ticks are from 1 to 65 inclusive.
intervals = unique([CLim(1) ceil(CLim(1)/Div)*Div : Div : floor(CLim(2)/Div)*Div CLim(2)]);
yt = interp1([CLim(1) CLim(2)], [1 65], intervals);
