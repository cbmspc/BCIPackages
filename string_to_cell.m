function c = string_to_cell(s,d);
% Convert a delimited string to a cell array
% E.g., input is    "blah 1" "blah 2", delimiter is ",  
%           output:    {'blah 1', 'blah 2'}

% Copyright 2001-2004 The MathWorks, Inc.
% $Revision: 1.1.8.1 $ $Date: 2004/07/21 06:23:56 $

c = {};
while containsValidString(s),
    [s1 s] = strtok(s, d);
    if containsValidString(s1)
        c = {c{:} s1};
    end
end

% ---------------------------------------
function ok = containsValidString(s)
% Decide whether there is still valid data in s.
% I.e., if s only contains separators, quotes, spaces,
% newlines, etc (in any combination), then it
% is not valid.
% This is to be decided in the context of 
% valid filenames, valid code symbols, etc.

goodChars = [ ...
    'abcdefghijklmnopqrstuvwxyz' ...
    'ABCDEFGHIJKLMNOPQRSTUVWXYZ' ...
    '1234567890' ...
    '_~-.!#$%'];
% !"#$%&'()*+,-./0123456789:;<=>?@
% [\]^_`
s2 = strtok(s, goodChars);
% If s2 does not contain any of these characters,
% s and s2 will be equal.
ok = ~isequal(s2, s);
