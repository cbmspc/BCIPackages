tmp_loadwork_indexfile = [getworkdir() filesep 'matlab-' crc32(lower(getusername)) '-' crc32(lower(pwd)) '-index' '.mat'];
tmp_loadwork_smallfile = [getworkdir() filesep 'matlab-' crc32(lower(getusername)) '-' crc32(lower(pwd)) '-small' '.mat'];
if exist(tmp_loadwork_indexfile, 'file')
    fprintf('Loading your workspace from %s :\n', getworkdir());
    pause(0.01);
    fprintf('  *Loading small variables..\n');
    if exist(tmp_loadwork_smallfile, 'file')
        load(tmp_loadwork_smallfile);
    end
    clear tmp_savework_filelist
    load(tmp_loadwork_indexfile, 'tmp_savework_filelist');
    if exist('tmp_savework_filelist','var')
        for tmp_loadwork_index = 1:size(tmp_savework_filelist,1)
            tmp_loadwork_filename = tmp_savework_filelist{tmp_loadwork_index,1};
            if exist(tmp_loadwork_filename, 'file')
                tmp_dirinfo = dir(tmp_loadwork_filename);
                if ~isempty(tmp_dirinfo) && strcmpi(sha1sum([tmp_dirinfo(1).datenum tmp_dirinfo(1).bytes]),tmp_savework_filelist{tmp_loadwork_index,2})
                    tmp_loadwork_filevars = whos('-file', tmp_savework_filelist{tmp_loadwork_index,1});
                    if ~isempty(tmp_loadwork_filevars)
                        fprintf('  *Loading variable %s (%.1f MiB)..\n', tmp_loadwork_filevars(1).name, tmp_loadwork_filevars(1).bytes/1024^2);
                    else
                        fprintf('  *Loading variable (%.1f MiB saved)..\n', tmp_loadwork_filevars(1).bytes/1024^2);
                    end
                    load(tmp_savework_filelist{tmp_loadwork_index,1});
                end
            end
        end
    end
    clear tmp_loadwork_*
    fprintf('Done.\n');
else
    fprintf('Nothing to load. There is no saved index file specific to the folder you are in.\n');
    clear tmp_loadwork_*
end