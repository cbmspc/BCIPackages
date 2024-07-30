function [x2, scalefactor, mu2, mu] = equalize_traces (x1, x2)
% Translate and scale x2 to match the levels of x1

mu = median(x1);
md = mad(x1,1);

mu2 = median(x2);
md2 = mad(x2,1);

scalefactor = md/md2;

x2 = (x2-mu2) * scalefactor + mu;
