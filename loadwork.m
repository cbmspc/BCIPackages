function loadwork()
if ~exist([getusername() '-workspace.mat'], 'file')
    fprintf('Cannot load from saved workspace because there is no saved file (%s) in the current folder:\n %s\n Suggestion: Change folder first.\n\n', [getusername() '-workspace.mat'], cd);
    return
end
if ~isempty(evalin('base','whos'))
    optsa.Interpreter = 'tex';
    optsa.WindowStyle = 'modal';
    optsa.Default='Merge variables';
    answer = questdlg('You already have some variables in the workspace. Loading from saved workspace could overwrite some of them.', 'Loading from saved workspace', 'Clear then Load', 'Merge variables', 'Don''t Load', optsa);
    switch answer
        case 'Clear then Load'
            evalin('base', 'clear');
            evalin('base', 'load([getusername() ''-workspace.mat'']);');
        case 'Merge variables'
            evalin('base', 'load([getusername() ''-workspace.mat'']);');
    end

else
    evalin('base', 'load([getusername() ''-workspace.mat'']);');
end