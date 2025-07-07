function loadfileprogressbar(fraction, currentFile, windowTitle)
persistent wb wbstarttime wbfirstline
mem = memory;
memtxt = ['Memory: ' addByteString(mem.MemUsedMATLAB) ' used, ' addByteString(mem.MaxPossibleArrayBytes) ' free'];
etatxt = 'ETA: estimating...';
if exist('currentFile','var') && (ischar(currentFile) || isstring(currentFile))
    wbfirstline = currentFile;
end
if fraction <= 0
    % Make a new one
    if ~isempty(wb) && ishandle(wb)
        delete(wb);
    end
    if ~exist('windowTitle','var') || (~ischar(windowTitle) && ~isstring(windowTitle))
        windowTitle = 'Loading...';
    end
    wb = waitbar(0, [wbfirstline 10 etatxt 10 '(' memtxt ')'], 'Name', windowTitle);
    set(findobj(wb,'Interpreter','tex'),'Interpreter','none');
    wbstarttime = tic;
    return
elseif fraction >= 1
    if ~isempty(wb) && ishandle(wb)
        delete(wb);
    end
    return
end

if isempty(wb) || ~ishandle(wb)
    %warning('Waitbar handle has not been created or has been closed or cleared. Remake one with 0 progress.');
    return
end

elapsedtime = toc(wbstarttime);
elapsedfraction = fraction;
remainingfraction = 1 - fraction;
remainingtime = elapsedtime / elapsedfraction * remainingfraction;
etatxt = ['ETA: ' addDuration(remainingtime)];
waitbar(fraction, wb, [wbfirstline 10 etatxt 10 '(' memtxt ')']);
