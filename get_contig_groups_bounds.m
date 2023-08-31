% Get contiguous groups as bounds
% 

function Bounds = get_contig_groups_bounds (x)
if islogical(x)
    x = find(x);
end
CGS = get_contig_groups_string(sort(x),',');
CG = string_to_cell(CGS,',');

Bounds = zeros(length(CG),2);

for i = 1:length(CG)
    b = regexp(CG{i}, '^(\d+):(\d+)$', 'tokens', 'once');
    if ~isempty(b)
        Bounds(i,:) = [str2double(b{1}) str2double(b{2})];
        continue
    end
    b = regexp(CG{i}, '^(\d+)$', 'tokens', 'once');
    if ~isempty(b)
        Bounds(i,:) = [str2double(b{1}) str2double(b{1})];
        continue
    end
    
end