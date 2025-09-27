function savework(username, clearaftersave, foreground)
if ~exist('username', 'var') || isempty(username) || ~ischar(username)
    username = getusername();
end
if ~exist('clearaftersave','var') || ~isscalar(clearaftersave)
    clearaftersave = false;
end
if ~exist('foreground','var') || ~isscalar(foreground)
    foreground = false;
end
username = sanitizefilename(username);
usernameifdifferent = username;
if nargin == 0
    usernameifdifferent = '';
end
if isfile([username '-workspace.mat'])
    d = dir([username '-workspace.mat']);
    if isscalar(d) && now - d.datenum < 5/86400 %#ok<*TNOW1>
        error('Cannot save now. There is another instance of MATLAB currently saving into the same file.');
    end
end
w = evalin('base','whos');
if isempty(w)
    fprintf('Workspace is empty. No need to save.\n');
    return
end
wsize = sum([w.bytes]);
if wsize > 20e9
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
fprintf('\n\n********** SAVING YOUR WORKSPACE **********\nThe time now is: %s. Saving workspace variables (~%s MB) to disk under your username %s\n', datetime, addThousandsCommaSeparators(ceil(wsize/1024^2*1.1)), username);
persistent etimefactor 
if isempty(etimefactor)
    etimefactor = 1/17187500/60/4;
end
etimeminutes = wsize*etimefactor;
if wsize > 300e6
    bgnote = 'Background transfer may still continue afterwards.';
    if foreground
        bgnote = '';
    end
    fprintf('Estimated %.0f minutes before you can resume working. %s\nSubsequent saves are also faster, because unchanged variables are not saved again.\n           PLEASE WAIT...', etimeminutes, bgnote);
end

ccd = cd;
t_start = tic;
if clearaftersave
    if foreground
        evalin('base', ['saveparts([cd filesep ' '''' username '-workspace.mat''],''-verbose'',''-foreground''); clear;']);
    else
        evalin('base', ['saveparts([cd filesep ' '''' username '-workspace.mat''],''-verbose''); clear;']);
    end
else
    if foreground
        evalin('base', ['saveparts([cd filesep ' '''' username '-workspace.mat''],''-verbose'',''-foreground'');']);
    else
        evalin('base', ['saveparts([cd filesep ' '''' username '-workspace.mat''],''-verbose'');']);
    end
end
etimeminutes = toc(t_start)/60;
etimefactor = etimeminutes/wsize;


donestat = 'SAVED';
if clearaftersave
    donestat = 'SAVED AND CLEARED';
end
fgwarn = ['✅Foreground process is done. Background transfer may still be in progress in a command prompt window.\n' ...
    '                             Do not disconnect or turn off computer until it is done.\n'''];
if foreground
    fgwarn = '';
end
fprintf(['\n' ...
    fgwarn...
    '📁Location: <a href="matlab:winopen(''%s'');">Open in File Explorer</a> or <a href="matlab:cd(''%s'');">Change folder in MATLAB.</a>\n' ...
    '  To load this, go to the same folder and type loadwork %s\n' ...
    '********** WORKSPACE %s **********\n\n'], ...
    ccd, ccd, usernameifdifferent, donestat);

return


