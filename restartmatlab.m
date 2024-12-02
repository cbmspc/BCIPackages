function restartmatlab(arg1)
fprintf('Restarting MATLAB...\n');
% matlabrc.m will handle the autorestore
try
    tmp_aa = pathdef;
    tmp_bb = regexp(tmp_aa, '([^;]+(?i)MATLAB(?-i)[^;]*[\\/]toolbox[^;]+);', 'tokens', 'once');
    tmp_cc = regexp(tmp_bb{1}, '^(.+)\\toolbox', 'tokens', 'once');
    tmp_mm = [tmp_cc{1} '\bin\matlab.exe'];
    if exist(tmp_mm, 'file')
        tmp_batch = [feature('logdir') filesep 'restartmatlab' num2str(feature('getpid')) '.bat'];
        tmp_fh = fopen(tmp_batch, 'w');
        fprintf(tmp_fh, '@ECHO OFF\n');
        fprintf(tmp_fh, 'ECHO A new instance of MATLAB will be launched when the current instance has\n');
        fprintf(tmp_fh, 'ECHO finished exiting. Please wait about 10 more seconds after this window\n');
        fprintf(tmp_fh, 'ECHO disappears on its own.\n');
        fprintf(tmp_fh, ':LOOP\n');
        fprintf(tmp_fh, '%s%s%s%stasklist /fi "PID eq %i" 2>NUL | %s%s%s%sfind /i /n "MATLAB.exe" >NUL\n', getenv('SYSTEMROOT'), filesep, 'system32', filesep, feature('getpid'), getenv('SYSTEMROOT'), filesep, 'system32', filesep);
        fprintf(tmp_fh, 'IF %%ERRORLEVEL%% == 0 GOTO LOOP\n');
        fprintf(tmp_fh, 'ECHO Launching MATLAB in the background...\n');
        if nargin == 1
            fprintf(tmp_fh, 'START "" /min "%s" -nosplash -r "%s"\n', tmp_mm, arg1);
        else
            fprintf(tmp_fh, 'START "" /min "%s" -nosplash\n', tmp_mm);
        end
        fprintf(tmp_fh, 'START "" /min cmd /c DEL "%s"\n', tmp_batch);
        fclose(tmp_fh);
        system(sprintf('start "" /D "%s" "cmd /c %s"', feature('logdir'), tmp_batch));
        clear tmp_batch tmp_fh tmp_aa tmp_bb tmp_cc tmp_mm
    end
catch ex
    fprintf('Auto-restart failed. MATLAB is exiting in about 10 seconds. Please restart manually.\n');
    ex.message
    pause(10);
end

try
    disp(getString(message('MATLAB:finishsav:SavingWorkspaceData')));
    evalin('base', 'save([feature(''logdir'') filesep ''matlabexitsave'' num2str(java.lang.System.currentTimeMillis) ''.mat''], ''-v7.3'', ''-nocompression'');');
    autosave delete
catch
end
quit force
