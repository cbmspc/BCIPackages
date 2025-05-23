% Standard error of the mean
function se = sem(x)
se = std(x,'omitmissing')/nnz(isfinite(x));

