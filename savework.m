tmp_savework_whos = whos; 
tmp_savework_bytes = cat(1,tmp_savework_whos.bytes);
[tmp_savework_bytes, tmp_savework_order] = sort(tmp_savework_bytes, 'descend');
fprintf('Saving your workspace (%.1f MiB) to %s :\n', sum(tmp_savework_bytes)/1024^2, getworkdir());
tmp_savework_filelist = cell(length(tmp_savework_bytes),2);
for tmp_savework_index = 1:length(tmp_savework_bytes)
    if tmp_savework_bytes(tmp_savework_index) >= 8*1024^2
        fprintf('  *Variable %s (%.1f MiB): ', tmp_savework_whos(tmp_savework_order(tmp_savework_index)).name, tmp_savework_whos(tmp_savework_order(tmp_savework_index)).bytes/1024^2);
        fprintf('Hashing.. ');
        tmp_savework_hash = sha1sum(eval(tmp_savework_whos(tmp_savework_order(tmp_savework_index)).name));
        tmp_savework_filename = [getworkdir() filesep 'matlab-' crc32(lower(getusername)) '-' crc32(lower(pwd)) '-' 'large-' crc32(lower(tmp_savework_whos(tmp_savework_order(tmp_savework_index)).name)) '-' tmp_savework_hash '.mat'];
        tmp_savework_filelist{tmp_savework_index,1} = tmp_savework_filename;
        if ~exist(tmp_savework_filename, 'file')
            fprintf('Saving.. ');
            save(tmp_savework_filename, tmp_savework_whos(tmp_savework_order(tmp_savework_index)).name, '-v7.3');
            fprintf('Saved.\n');
        else
            fprintf('Already saved.\n');
        end
    else
        break;
    end
end
pause(0.01);
fprintf('  *Saving small variables..\n');
tmp_savework_filename = [getworkdir() filesep 'matlab-' crc32(lower(getusername)) '-' crc32(lower(pwd)) '-' 'small' '.mat'];
save(tmp_savework_filename, tmp_savework_whos(tmp_savework_order(tmp_savework_index:end)).name, '-v7.3');
fprintf('  *Saving index..\n');
for tmp_savework_index = 1:size(tmp_savework_filelist,1)
    if ~isempty(tmp_savework_filelist{tmp_savework_index,1}) && exist(tmp_savework_filelist{tmp_savework_index,1}, 'file')
        tmp_dirinfo = dir(tmp_savework_filelist{tmp_savework_index,1});
        if ~isempty(tmp_dirinfo)
            tmp_savework_filelist{tmp_savework_index,2} = sha1sum([tmp_dirinfo(1).datenum tmp_dirinfo(1).bytes]);
        end
    end
end
tmp_savework_filename = [getworkdir() filesep 'matlab-' crc32(lower(getusername)) '-' crc32(lower(pwd)) '-' 'index' '.mat'];
save(tmp_savework_filename, 'tmp_savework_filelist', '-v7.3');
clear tmp_savework_*
fprintf('Done.\n');

