% Calculates the nearest-neighbor distance for every point

function [avgdist, nndist] = dataproc_func_nndist (X)

dd = squareform((pdist(X)));
dd(find(dd==0)) = Inf;
nndist = min(dd);
avgdist = median(nndist);


% N = size(X,1);
% for i = 1:N
%     dist = [];
%     for j = [1:i-1,i+1:N]
%         dist = [dist sqrt(sum((X(i,:) - X(j,:)).^2))];
%     end
%     lowestdist(i) = min(dist);
% end
% lowestdist = sort(lowestdist);
