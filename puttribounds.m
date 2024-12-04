% Opposite of getbounds
function Signal = puttribounds (TriBounds, SigLen)
Signal = putbounds(TriBounds, SigLen);

% if size(TriBounds,2) ~= 3
%     error('Expect input to be size N by 3, where each row is [startindex endindex numericlabel]');
% end
% 
% if ~exist('SigLen','var') || isempty(SigLen)
%     SigLen = max(TriBounds(:,2));
% end
% Signal = zeros(SigLen,1);
% 
% for i = 1:size(TriBounds,1)
%     Signal(TriBounds(i,1):TriBounds(i,2)) = TriBounds(i,3);
% end
% 
