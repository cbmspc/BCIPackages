% Robust estimate of standard deviation using median absolute deviation 
function sd = robuststd(x)
sd = mad(x,1)*1.4826;
