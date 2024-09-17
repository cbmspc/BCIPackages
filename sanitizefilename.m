function outputname = sanitizefilename (inputname)
% Remove illegal chars, space at beginning or end
outputname = regexprep(inputname, '^\s+|\s+$|[\x00-\x1F]', '');

% Replace illegal chars for filename
outputname = regexprep(outputname, '/|\\|?|%|*|:|\||"|<|>', '=');

% Replace space with + sign (also de-duplicates them)
outputname = regexprep(outputname, ' {1,}', '+');
