function savework(username, clearaftersave)
if ~exist('username', 'var') || isempty(username) || ~ischar(username)
    username = getusername();
end
if ~exist('clearaftersave','var') || ~isscalar(clearaftersave)
    clearaftersave = false;
end
username = sanitizefilename(username);
if isfile([username '-workspace.mat'])
    d = dir([username '-workspace.mat']);
    if isscalar(d) && now - d.datenum < 3/86400 %#ok<*TNOW1>
        error('Cannot save now. There is another instance of MATLAB currently saving into the same file.');
    end
end
w = evalin('base','whos');
if isempty(w)
    fprintf('Workspace is empty. No need to save.\n');
    return
end
wsize = sum([w.bytes]);
if wsize > 10e9
    optsa.Interpreter = 'tex';
    optsa.WindowStyle = 'modal';
    optsa.Default='Save them all';
    answer = questdlg([sprintf('Your workspace contains variables taking up ~%.1f GB space. Saving all of these variables can take ~%.1f minutes. Are you sure you want to save?', wsize/1024^3*1.1, wsize/10e6/60*1.1)], 'Saving the workspace', 'Save them all', 'Never mind', optsa);
    switch answer
        case 'Save them all'
        case 'Never mind'
            fprintf('You decided not to save the workspace.\n');
            return
        otherwise
            fprintf('You decided not to save the workspace.\n');
            return
    end

end
fprintf('\n\n********** SAVING YOUR WORKSPACE **********\nSaving workspace variables to disk. Please wait.\n');
%etimeminutes = wsize/10e6/60*1.1;
etimeminutes = wsize/17187500/60;
if wsize > 300e6
    fprintf('Because you are saving ~%s MB of data, this can take a while\n (estimated %.1f minutes, but some of this time is in the background where you can continue working)...', addThousandsCommaSeparators(ceil(wsize/1024^2*1.1)), etimeminutes);
    % try
    %     system(['explorer.exe "' cd '"']);
    % catch
    % end
end

ccd = cd;
if clearaftersave
    evalin('base', ['saveparts([cd filesep ' '''' username '-workspace.mat''],''-verbose''); clear;']);
else
    evalin('base', ['saveparts([cd filesep ' '''' username '-workspace.mat''],''-verbose'');']);
end

donestat = 'SAVED';
if clearaftersave
    donestat = 'SAVED AND CLEARED';
end
fprintf('\n✅Foreground process is done. (Background transfer may still be in progress in a command prompt window.)\nLocation: %s\nTo load this, go to the same folder and type loadwork\n********** WORKSPACE %s **********\n\n', ccd, donestat);

return


