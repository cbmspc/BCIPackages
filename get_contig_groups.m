function [ContigGroups, CLength] = get_contig_groups (x)
% function [ContigGroups, CLength] = get_contig_groups (x)
%
% Get the contiguous groups of integers from x as a cell array
% For example, if x = [1 2 3 4 9 10 11 12 20 21 100]
% Then ContigGroups = {[1 2 3 4], [9 10 11 12], [20 21], [100]}
%
% CLength returns the length of each group
%

if islogical(x)
    x = find(x);
end

x = [x(:).' inf];
ContigGroups = cell(1,length(x));
CLength = zeros(1,length(x));
k = 0;
while ~isempty(x)
    k = k+1;
    ContigGroups{k} = x(1:(find(diff(x) > 1 | diff(x) < 1)));
    CLength(k) = length(ContigGroups{k});
    x = x(1+(find(diff(x) > 1 | diff(x) < 1)):end);
end
ContigGroups = ContigGroups(1:k-1);
CLength = CLength(1:k-1);
