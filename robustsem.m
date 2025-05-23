% Robust estimate of standard error using median absolute deviation 
function se = robustsem(x)
se = robuststd(x)/nnz(isfinite(x));

