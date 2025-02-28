function saveworkclear()
if exist([getusername() '-workspace.mat'], 'file')
    d = dir([getusername() '-workspace.mat']);
    if now - d.datenum < 3/86400 %#ok<*TNOW1>
        error('Cannot save now. There is another instance of MATLAB currently saving into the same file.');
    end
end
w = evalin('base','whos');
if isempty(w)
    fprintf('Workspace is empty. No need to save.\n');
    return
end
wsize = sum([w.bytes]);
if wsize > 1e9
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
fprintf('\n\n********** SAVING YOUR WORKSPACE AND CLEARING WHEN SAVED **********\nSaving workspace variables to disk. Please wait.\n');
etimeminutes = wsize/10e6/60*1.1;
if wsize > 300e6
    fprintf('Because you are saving ~%s MB of data, this can take a while (estimated %.1f minutes)...', addComma(ceil(wsize/1024^2*1.1)), etimeminutes);
    % try
    %     system(['explorer.exe "' cd '"']);
    % catch
    % end
end

ccd = cd;
evalin('base', 'save([getusername() ''-workspace.mat''], ''-v7.3'', ''-nocompression''); clear;');

fprintf('\n Done! Saved in %s\n********** WORKSPACE SAVED AND CLEARED **********\n\n', ccd);

return


function numOut = addComma(numIn)
%https://www.mathworks.com/matlabcentral/answers/96131-is-there-a-format-in-matlab-to-display-numbers-such-that-commas-are-automatically-inserted-into-the
jf=java.text.DecimalFormat; % comma for thousands, three decimal places
numOut= char(jf.format(numIn)); % omit "char" if you want a string out
