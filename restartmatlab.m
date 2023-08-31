function restartmatlab(arg1)
if exist([gettmpdir() filesep 'activeworkspace.mat'],'file')
    tmp_a = dir([gettmpdir() filesep 'activeworkspace.mat']);
    if (now - tmp_a.datenum)*86400 < 30
        fprintf('An instance of MATLAB is being restarted.\nPlease wait up to 30 seconds and try again.\n');
        return
    end
end
fprintf('Restarting MATLAB after saving workspace and path...\n');
evalin('base','tmp_cwd = pwd;');
evalin('base','tmp_custompath = path();');
evalin('base','tmp_matlabtitle = getmatlabtitle();');
evalin('base','save([gettmpdir() filesep ''activeworkspace.mat'']);');
try
    a = pathdef;
    b = regexp(a, '([^;]+(?i)MATLAB(?-i)[^;]*[\\/]toolbox[^;]+);', 'tokens', 'once');
    c = regexp(b{1}, '^(.+)\\toolbox', 'tokens', 'once');
    m = [c{1} '\bin\matlab.exe'];
    if exist(m, 'file')
        cd(getbciprogramsdir());
        if nargin == 1
            eval(['!' '"' m '" -nosplash -r "' arg1 '"']);
        else
            eval(['!' '"' m '" -nosplash']);
        end
    end
catch
    fprintf('Auto-restart failed. MATLAB is exiting.\n');
    pause(5.0);
end
quit('force');
