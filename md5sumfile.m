function hash = md5sumfile (File)
[s, msg] = system(['md5sum -b "' File '"']);
if s == 0
    msg = regexp(msg, '([0-9a-z]{32}) ', 'tokens', 'once');
    hash = msg{1};
else
    hash = repmat('0', 1, 32);
end