% Convert signal to common average reference
% Channels with no deviation are ignored (e.g. disconnected channels)
function [out, m] = car (in)

if iscell(in)
    for i = length(in):-1:1
        [out{i}, m{i}] = car(in{i});
    end
    return
end


if size(in,3) > 1
    % 3D matrix (ch x time x trial)
    out = in;
    s = std(in,[],2);
    for tr = 1:size(in,3)
        m = mean(in(s(:,:,tr)>0,:,tr),1);
        for ch = find(s(:,:,tr).' > 0)
            out(ch,:,tr) = in(ch,:,tr) - m;
        end
    end
else
    % 2D matrix (time x ch)
    out = in;
    if size(in,2) > size(in,1)
        warning(['Input data may need to be transposed first (# of time points = ' num2str(size(in,1)) ', # of channels = ' num2str(size(in,2)) ').']);
    end
    s = std(in,[],1);
    m = mean(in(:,s>0),2);
    for ch = find(s > 0)
        out(:,ch) = in(:,ch) - m;
    end
end
