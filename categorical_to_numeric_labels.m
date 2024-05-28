function [n, u] = categorical_to_numeric_labels (labels)
c = categorical(labels);
u = categories(c);
n = nan(size(labels));
for i = 1:length(u)
    n(c == u{i}) = i;
end