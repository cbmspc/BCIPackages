function output = cellequal (cell1, cell2)
s1 = size(cell1);
s2 = size(cell2);

if numel(s1) ~= numel(s2)
    output = 0;
    return
end

if min(s1 == s2) == 0
    output = 0;
    return
end

for i = 1:numel(cell1)
    if ~isequal(cell1{i},cell2{i})
        output = 0;
        return
    end
end

output = 1;
