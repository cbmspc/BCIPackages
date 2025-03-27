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


parts_bytes_threshold = 32768;

if length(filename) < 2
    error('Filename is too short');
end
smallfilename = [filename fileext];
largefilenameprefix = regexprep([filename fileext], '\.mat$', '--');
largefilenamesuffix = '.mat';

variablenames = varargin;


whos_list = evalin(workspace,'whos');
if isempty(variablenames)
    whos_keeplist = whos_list;
else
    [~, ia] = intersect({whos_list.name},variablenames);
    whos_keeplist = whos_list(ia);
end

keeplist = {};

whos_smalllist = whos_keeplist([whos_keeplist.bytes] <= parts_bytes_threshold);
evalin(workspace, ['save(''' [pathname filesep smallfilename] ''',''' cell_to_string([{whos_smalllist.name} {'-mat'} {'-v7.3'} {'-nocompression'}],''',''') ''');']);
keeplist = [keeplist; {smallfilename}];

whos_largelist = whos_keeplist([whos_keeplist.bytes] > parts_bytes_threshold);
for i = 1:length(whos_largelist)
    whos_largelist(i).datahash = datahash(evalin(workspace,whos_largelist(i).name), struct('Input','array','Method','SHA-256','Format','hex'));
end

for i = 1:length(whos_largelist)
    dlist = dir([pathname filesep largefilenameprefix whos_largelist(i).datahash '-*' largefilenamesuffix]);
    % Check for tampering
    comparison_filename = [largefilenameprefix whos_largelist(i).datahash '-' datahash([dlist.bytes dlist.datenum],struct('Input','array','Method','MD5','Format','hex')) largefilenamesuffix];
    if isscalar(dlist) && strcmp(dlist.name, comparison_filename)
        keeplist = [keeplist; {dlist.name}]; %#ok<*AGROW>
        continue
    elseif length(dlist) > 1
        for j = 1:length(dlist)
            delete([dlist(j).folder filesep dlist(j).name]);
        end
    end
    savestr = sprintf('save(''%s'',''%s'',''-mat'',''-v7.3'',''-nocompression'');', [pathname filesep largefilenameprefix whos_largelist(i).datahash '-still_saving_please_wait' largefilenamesuffix], whos_largelist(i).name);
    evalin(workspace, savestr);
    dlist = dir([pathname filesep largefilenameprefix whos_largelist(i).datahash '-still_saving_please_wait' largefilenamesuffix]);
    comparison_filename = [largefilenameprefix whos_largelist(i).datahash '-' datahash([dlist.bytes dlist.datenum],struct('Input','array','Method','MD5','Format','hex')) largefilenamesuffix];
    movefile([pathname filesep largefilenameprefix whos_largelist(i).datahash '-still_saving_please_wait' largefilenamesuffix], [pathname filesep comparison_filename], 'f');
    keeplist = [keeplist; {comparison_filename}];
end

dlist = dir([pathname filesep largefilenameprefix '*' largefilenamesuffix]);
[~, ia] = setdiff({dlist.name}, keeplist);
dlist = dlist(ia);

for i = 1:length(dlist)
    delete([dlist(i).folder filesep dlist(i).name]);
end


