function [dpath, amask, atype] = getcccdatadir (varargin)

% options:
% -inlab: Assumes that your current location is within lab LAN
% -outlab: Assumes that your current location is outside lab LAN
% -offcampus: Assumes that you are off campus, and cannot connect to SMB
% -cache: Use the featureserver cache
% -mkdir: If the last child dir does not exist, automatically make it

vgin = varargin;
opts = {};
for i = 1:length(vgin)
    if length(vgin) > 1 && vgin{i}(1) == '-'
        if  vgin{i}(2) ~= '-'
            opts{i} = vgin{i}(2:end); %#ok<AGROW>
        else
            i = i + 1; %#ok<FXSET>
            break
        end
    else
        break
    end
end
vgin = vgin(i:end);

if ~isempty(vgin)
    subdir = cell_to_string(vgin, filesep);
else
    subdir = '';
end

[dpath, amask, atype] = getcccdatadir_main(subdir, opts);

if isempty(dpath) || strcmpi(subdir, 'managesubjects')
    if strcmpi(subdir, 'managesubjects')
        if exist('D:\managesubjects', 'dir')
            dpath = 'D:';
        elseif exist('D:\SnapMirror\managesubjects', 'dir')
            dpath= 'D:\SnapMirror';
        end
    elseif strcmpi(subdir, 'managesubjects\studies')
        if exist('D:\managesubjects\studies', 'dir')
            dpath = 'D:';
        elseif exist('D:\SnapMirror\managesubjects\studies', 'dir')
            dpath= 'D:\SnapMirror';
        end
    end
    
end


if exist('subdir', 'var') && ~isempty(subdir) && ischar(subdir)
    if exist([dpath filesep subdir], 'dir') || exist([dpath filesep subdir], 'file')
        dpath = [dpath filesep subdir];
    elseif ismember('mkdir', opts)
        subdir2 = regexprep(subdir, '[\\/][^\\/]+$', '');
        [dpath, amask, atype] = getcccdatadir_main(subdir2, opts);
        if exist([dpath filesep subdir2], 'dir') || exist([dpath filesep subdir2], 'file')
            dpath = [dpath filesep subdir];
            mkdir(dpath);
        else
            error('Specified directory %s and its parent %s do not exist.', subdir, subdir2);
        end
    else
        error('Specified directory %s does not exist.', subdir);
    end
end

function [dpath, amask, atype] = getcccdatadir_main (searchfile, opts)

global g_CCC_DATA_DIR_OVERRIDE
persistent directorypath accessmask accesstype searchedfile optedstr

dpath = '';
amask = 0;
atype = 0;

% accesstype can be 
% 3 for local and local-like access
% 2 for remote access
% 0 for no access

% accessmask can be
% 3 for unrestricted read-write access
% 2 for read-write access to files and dirs
% 1 for read only access to files and dirs
% 0 otherwise

pc{1,1} = {
    'D:'
    '\\ucibciserver.stemcell.uci.edu\Research_ncc135p'
    'D:\SnapMirror\ncc135p_data'
    'D:\SnapMirror\ncc11bw_data'
    'E:'
    'F:'
    };

pc{2,1} = {
    };

choices = true(2,1);

if ismember('cache', opts)
    choices(1) = 0;
end
if ismember('inlab', opts)
    choices(2) = 0;
elseif ismember('outlab', opts)
    choices(2) = 0;
elseif ismember('offcampus', opts)
    choices(2) = 0;
end

searchorder.pc = cat(1,pc{choices});

searchorder.unix = {
    '/home/.Backups/ccc_data/current'
};

if exist('g_CCC_DATA_DIR_OVERRIDE', 'var') && ~isempty(g_CCC_DATA_DIR_OVERRIDE) && exist(g_CCC_DATA_DIR_OVERRIDE, 'dir')
    searchorder.pc = {g_CCC_DATA_DIR_OVERRIDE};
    searchorder.unix = {g_CCC_DATA_DIR_OVERRIDE};
end


searchcheck.dir = {
    'experiment'
    %'managesubjects'
    %'bci_ProgramData'
    'intermediate'
    'rawoutput'
};

searchcheck.file = {
};

if exist('searchfile','var') && ~isempty(searchfile) && ischar(searchfile)
    searchcheck.searchfile = searchfile;
else
    searchfile = '';
end

% if isempty(computername)
%     computername = getcomputername();
% end

optstr = cell_to_string(opts, '-');

if ~isempty(directorypath) && ~isempty(accessmask) && ~isempty(accesstype) && ~isempty(searchedfile) && ~isempty(optedstr)
    found = validate(searchcheck, directorypath);
    if found && strcmp(searchfile, searchedfile) && strcmp(optstr, optedstr)
        dpath = directorypath;
        amask = accessmask;
        atype = accesstype;
        return
    end
end


if ispc
    for i = 1:length(searchorder.pc)
        if 1 || exist(searchorder.pc{i},'dir')
            [found, accessmask, accesstype] = validate(searchcheck, searchorder.pc{i});
            if found
                directorypath = searchorder.pc{i};
                dpath = directorypath;
                amask = accessmask;
                atype = accesstype;
                searchedfile = searchfile;
                optedstr = optstr;
                return
            end
        end
    end
else
    for i = 1:length(searchorder.unix)
        if exist(searchorder.unix{i},'dir')
            [found, accessmask, accesstype] = validate(searchcheck, searchorder.unix{i});
            if found
                directorypath = searchorder.unix{i};
                dpath = directorypath;
                amask = accessmask;
                atype = accesstype;
                searchedfile = searchfile;
                optedstr = optstr;
                return
            end
        end
    end
end


function [found, accessmask, accesstype] = validate (searchcheck, searchpath)
found = 1;
for j = 1:length(searchcheck.dir)
    if exist([searchpath filesep searchcheck.dir{j}],'dir')
        found = found*1;
    else
        found = found*0;
    end
end
for j = 1:length(searchcheck.file)
    if exist([searchpath filesep searchcheck.file{j}],'file')
        found = found*1;
    else
        found = found*0;
    end
end

if isfield(searchcheck, 'searchfile')
    if exist([searchpath filesep searchcheck.searchfile],'file') || exist([searchpath filesep searchcheck.searchfile],'dir')
        found = found*1;
    else
        found = found*0;
    end
end

if nargout >= 2
    accessmask = 0;
    a = fileaccesscheck(searchpath);
    if min(a(1:2))
        % has read access
        accessmask = 1;
        if min(a(3:7))
            % also has modify access
            accessmask = 2;
            if min(a(8:9))
                % also has modifyattribute access
                accessmask = 3;
            end
        end
    end
end

if nargout >= 3
    accesstype = 0;
    if length(searchpath) >= 2
        accesstype = 3;
        if strcmp(searchpath(1:2),'\\') || strcmp(searchpath(1:2),'//')
            accesstype = 2;
        end
    end
end
