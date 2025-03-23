function str = addSpacesBeforeString(str, totalLength)
L = length(str);
str = [repmat(' ',[1,totalLength-L]) str];
