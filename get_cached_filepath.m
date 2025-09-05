function filepath = get_cached_filepath(original_filepath, skipIfNotCached)
disk_space_cleanup_threshold = 512*1024^3;

original_filepath = lower(original_filepath);
alternate_filepath = original_filepath;

alternate_filepath = regexprep(alternate_filepath, '^\\\\ucibciserver\.d\.zind\.org\\', '\\\\ucibciserver\.stemcell\.uci\.edu\\');
alternate_filepath = regexprep(alternate_filepath, '^\\\\ucibciserver\.myqnapcloud\.com\\', '\\\\ucibciserver\.stemcell\.uci\.edu\\');
alternate_filepath = regexprep(alternate_filepath, '^\\\\128\.200\.248\.115\\', '\\\\ucibciserver\.stemcell\.uci\.edu\\');
alternate_filepath = regexprep(alternate_filepath, '^\\\\10\.141\.150\.188\\', '\\\\ucibciserver\.stemcell\.uci\.edu\\');

if strcmpi(getcomputername(), 'NCC18JQ')
    alternate_filepath = regexprep(alternate_filepath, '^\\\\ucibciserver\.stemcell\.uci\.edu\\', '\\\\10\.141\.150\.188\\');
end

if isfolder('D:\rosync') || isfolder('P:\rosync')
    alternate_filepath_rootfolder = regexp(alternate_filepath, '^\\\\ucibciserver\.stemcell\.uci\.edu\\([^\\]+)','tokens','once');
    if ~isempty(alternate_filepath_rootfolder)
        alternate_filepath_rootfolder = alternate_filepath_rootfolder{1};
        if isfolder(['D:\rosync\' alternate_filepath_rootfolder])
            alternate_filepath = regexprep(alternate_filepath, ['\\\\ucibciserver\.stemcell\.uci\.edu\\' alternate_filepath_rootfolder '\\'], ['D:\\rosync\\' alternate_filepath_rootfolder '\\']);
        elseif isfolder(['P:\rosync\' alternate_filepath_rootfolder])
            alternate_filepath = regexprep(alternate_filepath, ['\\\\ucibciserver\.stemcell\.uci\.edu\\' alternate_filepath_rootfolder '\\'], ['P:\\rosync\\' alternate_filepath_rootfolder '\\']);
        end
    end
end

dinfo = dir(original_filepath);
adinfo = dir(alternate_filepath);
if isempty(dinfo) || numel(dinfo) > 1 || dinfo.isdir
    % Not a file, so we don't care
    filepath = original_filepath;
    return
end
if isempty(adinfo) || numel(adinfo) > 1 || adinfo.isdir || adinfo.datenum ~= dinfo.datenum || adinfo.bytes ~= dinfo.bytes
    alternate_filepath = original_filepath;
end

% These original paths are equivalent (same IP address) and are converted
% to the canonical for hashing purpose
dinfo.folder = get_cached_filepath_getCanonicalFolderNames(dinfo.folder);



dotfile_extension = regexp(dinfo.name, '\.[^.]+$', 'match', 'once');

opt.Input = 'array';
opt.Format = 'hex';
opt.Method = 'SHA-256';

hashed_dinfo = datahash(dinfo, opt);
if length(hashed_dinfo) < 32
    % output of datahash is malformed (expect 32 characters for MD5 (please
    % dont use anything cheaper than MD5))
    filepath = original_filepath;
    return
end
hashed_dinfo_first4 = hashed_dinfo(1:4);

if ~exist('skipIfNotCached','var') || isempty(skipIfNotCached) || ~isscalar(skipIfNotCached)
    skipIfNotCached = false;
end

cachedir = get_nascache_dir();

if startsWith(original_filepath,lower(cachedir))
    % It is already a NAS cache
    filepath = original_filepath;
    return
end


cachesubdir = [cachedir filesep hashed_dinfo_first4];
if ~isfolder(cachesubdir)
    mkdir(cachesubdir);
end
%cachesubdir = cachedir;
cachefile = [cachesubdir filesep hashed_dinfo dotfile_extension];
lastaccessedfile = [cachesubdir filesep hashed_dinfo '-la.tmp'];

if ~isfolder(cachesubdir)
    % Can't create the cache dir. Something is wrong
    filepath = original_filepath;
    return
end

FileObj = java.io.File(cachesubdir);
usable_bytes = FileObj.getUsableSpace;
usable_bytes = usable_bytes - dinfo.bytes;
deficient_bytes = 0;
if usable_bytes < disk_space_cleanup_threshold
    % Clear more so we don't need to check again so soon
    deficient_bytes = floor(disk_space_cleanup_threshold*1.50 - usable_bytes); 
end

persistent lastchecked
if ~isscalar(lastchecked)
    lastchecked = 0;
end

if isfolder(cachedir) && deficient_bytes > 0 && now - lastchecked > 0.5 %#ok<*TNOW1>
    lastchecked = now;
    list_cached = dir([cachedir filesep '**' filesep '*']);
    list_cached = list_cached(~[list_cached.isdir]);
    list_cached = list_cached(cellfun(@isempty,regexp({list_cached.name}, '^\.\.?$', 'match', 'once')));
    list_files = list_cached(cellfun(@isempty,regexp({list_cached.name}, '-la\.tmp$', 'match', 'once')));
    list_cached = list_cached(~cellfun(@isempty,regexp({list_cached.name}, '-la\.tmp$', 'match', 'once')));
    ld = [list_cached.datenum];
    [~, ind] = sort(ld);
    list_cached = list_cached(ind);
    clear ind ld
    lb = [list_cached.bytes];
    lb_cum = cumsum(lb);
    if deficient_bytes > 0
        last_to_axe_by_size = find(lb_cum >= deficient_bytes,1);
        axethese = regexprep({list_cached(1:last_to_axe_by_size).name}, '-la.tmp$', '');
        if ~isempty(axethese)
            for i = 1:length(axethese)
                if length(axethese{i}) < 32
                    warning('Something important is broken!');
                    break
                end
                delete([cachesubdir filesep axethese{i} '.*']);
                delete([cachesubdir filesep axethese{i} '-la.tmp']);
            end
        end
    end

    % Delete any orphaned files (files without matching la files)
    list_files_la = regexprep({list_files.name}, '\..*$', '-la.tmp');
    for i = 1:length(list_files)
        if ~isfile([list_files(i).folder filesep list_files_la{i}])
            delete([list_files(i).folder filesep list_files(i).name]);
        end
    end
end


if isfolder(cachesubdir) && isfile(cachefile)
    filepath = cachefile;
    linfo = dir(cachefile);
    if linfo.bytes ~= dinfo.bytes || linfo.datenum ~= dinfo.datenum
        % Update immediately
        copystatus = copyfile(alternate_filepath, cachefile, 'f');
        if ~copystatus
            copystatus = copyfile(original_filepath, cachefile, 'f');
        end
        if ~copystatus
            % Can't copy the file. Something is wrong
            filepath = original_filepath;
            return
        end
    end
elseif skipIfNotCached
    filepath = original_filepath;
else
    % We cache it right now!
    if ~isfolder(cachesubdir)
        mkdir(cachesubdir);
    end
    if ~isfolder(cachesubdir)
        % Can't create the cache dir. Something is wrong
        filepath = original_filepath;
        return
    end
    copystatus = copyfile(alternate_filepath, cachefile, 'f');
    if ~copystatus
        copystatus = copyfile(original_filepath, cachefile, 'f');
    end
    if copystatus
        filepath = cachefile;
    else
        % Can't copy the file. Something is wrong
        filepath = original_filepath;
    end
    
end

% Update the last-accessed metafile
try fclose(fopen(lastaccessedfile, 'w')); end %#ok<TRYNC>

