function numOut = addThousandsCommaSeparators(numIn)
%https://www.mathworks.com/matlabcentral/answers/96131-is-there-a-format-in-matlab-to-display-numbers-such-that-commas-are-automatically-inserted-into-the
jf=java.text.DecimalFormat; % comma for thousands, three decimal places
numOut= char(jf.format(numIn)); % omit "char" if you want a string out
