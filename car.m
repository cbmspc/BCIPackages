% Convert signal to common average reference
% Channels with no deviation are ignored (e.g. disconnected channels)
function [out, m] = car (in)
out = in;
s = std(in,[],1);
m = mean(in(:,s>0),2);
for i = find(s > 0)
    out(:,i) = in(:,i) - m;
end
