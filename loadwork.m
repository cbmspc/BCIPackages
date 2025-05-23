function loadwork(username)
if ~exist('username', 'var') || isempty(username) || ~ischar(username)
    username = getusername();
end
username = sanitizefilename(username);
if isfile([username '-workspace.mat'])
    filename = [username '-workspace.mat'];
elseif isfile('matlab.mat')
    filename = 'matlab.mat';
else
    filename = '';
end
if isempty(filename)
    fprintf('Cannot load from saved workspace because there is no saved file in the current folder:\n %s\n Suggestion: Change folder first.\n\n', cd);
    return
end
if ~isempty(evalin('base','whos'))
    optsa.Interpreter = 'tex';
    optsa.WindowStyle = 'modal';
    optsa.Default='Merge variables';
    answer = questdlg('You already have some variables in the workspace. Loading from saved workspace could overwrite some of them.', 'Loading from saved workspace', 'Clear then Load', 'Merge variables', 'Don''t Load', optsa);
    switch answer
        case 'Clear then Load'
            fprintf('Clearing workspace and then loading...\n');
            evalin('base', 'clear');
            evalin('base', ['loadparts(''' cd filesep filename ''');']);
            fprintf('Loaded all variables from %s.\n', [cd filesep filename]);
        case 'Merge variables'
            fprintf('Loading...\n');
            evalin('base', ['loadparts(''' cd filesep filename ''');']);
            fprintf('Loaded all variables from %s.\n', [cd filesep filename]);
    end

else
    fprintf('Loading...\n');
    evalin('base', ['loadparts(''' cd filesep filename ''');']);
    fprintf('Loaded all variables from %s.\n', [cd filesep filename]);
end