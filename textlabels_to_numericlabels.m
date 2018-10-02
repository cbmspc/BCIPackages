% Convert text labels (cell array of string) to numeric labels 
% (number array and sorted label names)

function [numericlabels, labelnames] = textlabels_to_numericlabels (textlabels)
labelnames = unique(textlabels);
numericlabels = zeros(size(textlabels));
for c = 1:numel(labelnames)
    numericlabels(strcmp(textlabels, labelnames{c})) = c;
end
