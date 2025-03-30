function finish()
listwhos = evalin('base','whos');
if length(listwhos) < 2
    return
end
if sum([listwhos.bytes]) <= 104857600 && isfile(getappdata(0, 'g_matlabbaseworkspaceautosavefile'))
    % Up to 100 MB: Just autosave and quit.
    try 
        saveparts(getappdata(0, 'g_matlabbaseworkspaceautosavefile'), listwhos.name);
        return
    catch
    end
end
optsa.Interpreter = 'tex';
optsa.WindowStyle = 'modal';
optsa.Default='Don''t exit';
answer = questdlg(sprintf('\\fontsize{20}There are %i variables (%.2f MB) in the base workspace.\n\\color{black}Exit MATLAB \\color{orange}without \\color{black}saving?', length(listwhos), sum([listwhos.bytes])/1024^2), 'Exiting MATLAB', 'Exit without saving', 'Don''t exit', optsa);
switch answer
    case 'Exit without saving'
        return
    otherwise
        error('User decided not to exit MATLAB.');

end



