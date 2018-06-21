% function s = substr (str, sep, n)
% nth substring of str, separated by sep
function s = substr (str, sep, n)
s = string_to_cell(str, sep);
s = s{n};
