function filepath = get_cached_filepath(original_filepath, skipIfNotCached)
max_num_cached_files = 100000;
max_days_cached_files = 100;
disk_space_cleanup_threshold = 512*1024^3;

original_filepath = lower(original_filepath);
alternate_filepath = original_filepath;

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
%hashed_dinfo_first4 = hashed_dinfo(1:4);

if ~exist('skipIfNotCached','var') || isempty(skipIfNotCached) || ~isscalar(skipIfNotCached)
    skipIfNotCached = false;
end

cachedir = get_nascache_dir();
%cachesubdir = [cachedir filesep hashed_dinfo_first4];
cachesubdir = cachedir;
cachefile = [cachesubdir filesep hashed_dinfo dotfile_extension];
lastaccessedfile = [cachesubdir filesep hashed_dinfo '-la.tmp'];

if ~exist(cachedir, 'dir')
    % Can't create the cache dir. Something is wrong
    filepath = original_filepath;
    return
end

FileObj = java.io.File(cachesubdir);
usable_bytes = FileObj.getUsableSpace;
usable_bytes = usable_bytes - dinfo.bytes;
deficient_bytes = 0;
if usable_bytes < disk_space_cleanup_threshold
    deficient_bytes = disk_space_cleanup_threshold - usable_bytes;
end

if exist(cachesubdir, 'dir')
    list_cached = dir([cachesubdir filesep '*-la.tmp']);
    numfiles = length(list_cached);
    ld = [list_cached.datenum];
    [ld, ind] = sort(ld);
    list_cached = list_cached(ind);
    clear ind
    lb = [list_cached.bytes];
    lb_cum = cumsum(lb);
    if numfiles >= max_num_cached_files || any(now - ld > max_days_cached_files) || deficient_bytes > 0 %#ok<*TNOW1>
        last_to_axe_by_count = numfiles - max_num_cached_files + 1;
        last_to_axe_by_age = find(ld >= now - max_days_cached_files,1) - 1;
        last_to_axe_by_size = find(lb_cum >= deficient_bytes,1);
        numtoaxe = max([last_to_axe_by_count, last_to_axe_by_age, last_to_axe_by_size]);
        axethese = regexprep({list_cached(1:numtoaxe).name}, '-la.tmp$', '');
        if ~isempty(axethese)
            rstate = recycle;
            recycle('on');
            for i = 1:length(axethese)
                if length(axethese{i}) < 32
                    warning('Something important is broken!');
                    break
                end
                delete([cachesubdir filesep axethese{i} '.*']);
                delete([cachesubdir filesep axethese{i} '-la.tmp']);
            end
            recycle(rstate);
        end
    end
end


if exist(cachesubdir, 'dir') && exist(cachefile, 'file')
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
    if ~exist(cachesubdir, 'dir')
        mkdir(cachesubdir);
    end
    if ~exist(cachesubdir, 'dir')
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

