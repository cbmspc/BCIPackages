function s = cell_to_string (c, d)
s = '';
l = length(c);
if l == 0
    return
end
s = c{1};
for i = 2:l
    s = [s d c{i}];
end
