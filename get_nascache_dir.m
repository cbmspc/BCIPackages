function d = get_nascache_dir()
d = '';
common_cachedir = 'D:\ProgramData\AppData\Local\NASCache\CacheStorage';

if exist(common_cachedir, 'dir' )
    d = common_cachedir;
    return;
end

computername = getcomputername();
username = getusername();
d_users = ['D:\Users_' computername];
user_d_drive_folder = [d_users filesep username];
nascache = [user_d_drive_folder filesep 'AppData' filesep 'Local' filesep 'NASCache' filesep 'CacheStorage'];

if exist(user_d_drive_folder, 'dir')
    if ~exist(nascache, 'dir')
        mkdir(nascache);
    end
end
if exist(nascache, 'dir')
    d = nascache;
    return
end

