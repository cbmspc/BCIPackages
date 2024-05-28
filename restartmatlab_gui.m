function restartmatlab_gui(handles_to_disable, relaunchcommand)
options.Default = 'No, go back';
options.Interpreter = 'tex';
options.WindowStyle = 'modal';
relaunchprompt = '';
if nargin >= 1
    relaunchprompt = 'and relaunch this interface. ';
end
if nargin < 2
    relaunchcommand = '';
end
button = questdlg(['\color{red}\fontsize{13}\fontname{Arial Black}CONFIRM RESTART MATLAB \color{black}\fontsize{13}\fontname{Calibri}' 10 ...
    'MATLAB will close and restart itself after' 10 ...
    'about 15 seconds without showing anything' 10 ...
    relaunchprompt ...
    'Proceed?'], 'Confirm Restart MATLAB', 'Save workspace and restart', 'Clear workspace and restart', 'No, go back', options);
switch button
    case 'Save workspace and restart'
        if nargin >= 1
            try
                for i = 1:length(handles_to_disable)
                    set(handles_to_disable(i), 'Enable', 'off');
                end
                drawnow
            end
        end
        evalin('base',['restartmatlab(' '''' relaunchcommand '''' ')']);
    case 'Clear workspace and restart'
        if nargin >= 1
            try
                for i = 1:length(handles_to_disable)
                    set(handles_to_disable(i), 'Enable', 'off');
                end
                drawnow
            end
        end
        evalin('base', 'clear');
        evalin('base',['restartmatlab(' '''' relaunchcommand '''' ')']);
    otherwise
end
