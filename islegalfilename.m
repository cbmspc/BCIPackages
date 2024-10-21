function output = islegalfilename (filename)
pat = '[A-Za-z0-9!@#$%&()_+=,.;''~` -]';
output = isempty(regexprep(filename, pat, ''));

