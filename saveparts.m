% Save workspace variables to disk, but save large variables in separate
% files to speed up.
%
% filepath must be an absolute full path to a file ending in .mat

function saveparts(filepath, varargin)
workspace = 'base';
if ~endsWith(filepath, '.mat')
    error('Filepath must end with .mat');
end
[pathname, filename, fileext] = fileparts(filepath);
if ~isfolder(pathname)
    error('Not a folder: %s', pathname);
end
if ~strcmp(fileext,'.mat')
    error('Filepath must end with .mat');
end


parts_bytes_threshold = 1024^2;

if length(filename) < 2
    error('Filename is too short');
end
smallfilename = [filename fileext];
largefilenameprefix = regexprep([filename fileext], '\.mat$', '--');
largefilenamesuffix = '.mat';

variablenames = varargin;

verbose = false;
if ~isempty(variablenames) && iscell(variablenames)
    if strcmpi(variablenames{1},'-verbose')
        verbose = true;
    end
    variablenames = variablenames(cellfun(@isvarname,variablenames));
end


whos_list = evalin(workspace,'whos');
if isempty(variablenames)
    whos_keeplist = whos_list;
else
    [~, ia] = intersect({whos_list.name},variablenames);
    whos_keeplist = whos_list(ia);
end

keeplist = {};
if verbose
    wb = waitbar(0, '', 'Name', ['saveparts: ' smallfilename]);
    set(findobj(wb,'Interpreter','tex'),'Interpreter','none');
else
    wb = -1;
end
whos_smalllist = whos_keeplist([whos_keeplist.bytes] <= parts_bytes_threshold);
bytestracker = 0;
bytestrackermax = sum([whos_keeplist.bytes]);
if ~isempty(whos_smalllist)
    if ishandle(wb)
        waitbar(0, wb, sprintf('Saving the index file...'));
    end
    if sum([whos_smalllist.bytes]) < 2e31
        evalin(workspace, ['save(''' [pathname filesep smallfilename] ''',''' cell_to_string([{whos_smalllist.name} {'-mat'} {'-v7'}],''',''') ''');']);
    else
        evalin(workspace, ['save(''' [pathname filesep smallfilename] ''',''' cell_to_string([{whos_smalllist.name} {'-mat'} {'-v7.3'}],''',''') ''');']);
    end
else
    fclose(fopen([pathname filesep smallfilename],'w'));
end
keeplist = [keeplist; {smallfilename}];
bytestracker = bytestracker + sum([whos_smalllist.bytes]);

HashBytesPerSec = 55e6;
SaveBytesPerSec = 25e6;

builtin_basictypes = {'double', 'single', 'uint64', 'uint32', 'uint16', 'uint8', 'int64', 'int32', 'int16', 'int8', 'char', 'logical'};

whos_largelist = whos_keeplist([whos_keeplist.bytes] > parts_bytes_threshold);
for i = 1:length(whos_largelist)
    if ishandle(wb)
        waitbar(bytestracker/bytestrackermax, wb, sprintf('Hashing the variable: %s (%s MiB)\nTime to hash~ %s.', whos_largelist(i).name, addThousandsCommaSeparators(ceil(whos_largelist(i).bytes/1024^2)),estimate_eta(whos_largelist(i).bytes,HashBytesPerSec)));
    end
    bestinputtype = 'array';
    if whos_largelist(i).bytes < 2^31 && ~ismember(whos_largelist(i).class,builtin_basictypes)
        bestinputtype = 'tempfile';
    end
    tstart = tic;
    whos_largelist(i).datahash = [sanitizefilename(whos_largelist(i).name) '-' datahash(evalin(workspace,whos_largelist(i).name), struct('Input',bestinputtype,'Method','MD5','Format','hex'))];
    telap = toc(tstart);
    HashBytesPerSec = whos_largelist(i).bytes / telap;
    bytestracker = bytestracker + whos_largelist(i).bytes/3;
end

for i = 1:length(whos_largelist)
    dlist = dir([pathname filesep largefilenameprefix whos_largelist(i).datahash '-*' largefilenamesuffix]);
    % Check for tampering
    if ~isscalar(dlist) || isempty(dlist.datenum) || ~isnumeric(dlist.datenum)
        comparison_filename = 'not valid';
    else
        comparison_filename = [largefilenameprefix whos_largelist(i).datahash '-' datahash([dlist.bytes round(dlist.datenum*21600)],struct('Input','array','Method','MD5','Format','hex')) largefilenamesuffix];
    end
    if isscalar(dlist) && strcmp(dlist.name, comparison_filename)
        keeplist = [keeplist; {dlist.name}]; %#ok<*AGROW>
        continue
    elseif length(dlist) > 1
        for j = 1:length(dlist)
            delete([dlist(j).folder filesep dlist(j).name]);
        end
    end
    savestr = sprintf('save(''%s'',''%s'',''-mat'',''-v7.3'');', [pathname filesep largefilenameprefix whos_largelist(i).datahash '-still_saving_please_wait' largefilenamesuffix], whos_largelist(i).name);
    if whos_largelist(i).bytes < 2^31
        savestr = sprintf('save(''%s'',''%s'',''-mat'',''-v7'');', [pathname filesep largefilenameprefix whos_largelist(i).datahash '-still_saving_please_wait' largefilenamesuffix], whos_largelist(i).name);
    end
    if ishandle(wb)
        %waitbar(bytestracker/bytestrackermax, wb, sprintf('Saving the variable: %s (%s MiB)', whos_largelist(i).name, addThousandsCommaSeparators(ceil(whos_largelist(i).bytes/1024^2))));
        waitbar(bytestracker/bytestrackermax, wb, sprintf('Saving the variable: %s (%s MiB)\nTime to save~ %s.', whos_largelist(i).name, addThousandsCommaSeparators(ceil(whos_largelist(i).bytes/1024^2)),estimate_eta(whos_largelist(i).bytes,SaveBytesPerSec)));
    end
    tstart = tic;
    evalin(workspace, savestr);
    bytestracker = bytestracker + whos_largelist(i).bytes/3;
    dlist = dir([pathname filesep largefilenameprefix whos_largelist(i).datahash '-still_saving_please_wait' largefilenamesuffix]);
    comparison_filename = [largefilenameprefix whos_largelist(i).datahash '-' datahash([dlist.bytes round(dlist.datenum*21600)],struct('Input','array','Method','MD5','Format','hex')) largefilenamesuffix];
    movefile([pathname filesep largefilenameprefix whos_largelist(i).datahash '-still_saving_please_wait' largefilenamesuffix], [pathname filesep comparison_filename], 'f');
    keeplist = [keeplist; {comparison_filename}];
    telap = toc(tstart);
    SaveBytesPerSec = whos_largelist(i).bytes / telap;
    bytestracker = bytestracker + whos_largelist(i).bytes/3;
end

dlist = dir([pathname filesep largefilenameprefix '*' largefilenamesuffix]);
[~, ia] = setdiff({dlist.name}, keeplist);
dlist = dlist(ia);

for i = 1:length(dlist)
    if ishandle(wb)
        waitbar((i-1)/length(dlist), wb, sprintf('Cleaning up extraneous savefiles'));
    end
    delete([dlist(i).folder filesep dlist(i).name]);
end
if ishandle(wb)
    delete(wb);
end


function str = estimate_eta(bytes, BytesPerSec)
sec = bytes / BytesPerSec;
str = sprintf('%.0f min %.0f sec', floor(sec/60), mod(sec,60));

