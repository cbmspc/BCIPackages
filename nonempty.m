% Return only the non-empty elements of a cell array
function c = nonempty(c)
c = c(~ismissing(c));
