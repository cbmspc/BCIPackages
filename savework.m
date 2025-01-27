function savework()
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
    answer = questdlg([sprintf('Your workspace contains variables taking up %.1f GB space. Saving all of these variables can take %.1f minutes. Are you sure you want to save?', wsize/1024^3, wsize/10e6/60)], 'Saving the workspace', 'Save them all', 'Never mind', optsa);
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
fprintf('Saving workspace variables to disk... \n');
if wsize > 300e6
    fprintf('Because you are saving %s KB of data, this can take a while (estimated %.1f minutes). Please wait.\n You can watch the progress in the File Explorer (automatically opened for you).\n', addComma(ceil(wsize/1024)), wsize/10e6/60);
    try
        system(['explorer.exe "' cd '"']);
    catch
    end
end

evalin('base', 'save([getusername() ''-workspace.mat''], ''-v7.3'', ''-nocompression'');');

fprintf('Done!\n');

return


function numOut = addComma(numIn)
%https://www.mathworks.com/matlabcentral/answers/96131-is-there-a-format-in-matlab-to-display-numbers-such-that-commas-are-automatically-inserted-into-the
jf=java.text.DecimalFormat; % comma for thousands, three decimal places
numOut= char(jf.format(numIn)); % omit "char" if you want a string out
