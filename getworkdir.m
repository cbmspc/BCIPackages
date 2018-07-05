% This is currently hardcoded to D:\Users\YourUserName\work
function D = getworkdir ()
U = getusername();
S = filesep;
D = sprintf('D:%sUsers%s%s%swork', S, S, U, S);
if ~exist(D, 'dir')
    try
        mkdir(D);
    catch
        D = gettmpdir();
    end
end
