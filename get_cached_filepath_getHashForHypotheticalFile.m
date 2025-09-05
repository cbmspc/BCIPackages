function hash = get_cached_filepath_getHashForHypotheticalFile(fileName, folder, dateNum, bytes)
fileName = lower(fileName);
folder = lower(folder);
folder = get_cached_filepath_getCanonicalFolderNames(folder);
dinfo = struct('name', fileName, 'folder', folder, 'date', datestr(dateNum), 'bytes', bytes, 'isdir', false, 'datenum', dateNum); %#ok<DATST>

opt.Input = 'array';
opt.Format = 'hex';
opt.Method = 'SHA-256';
hash = datahash(dinfo, opt);
