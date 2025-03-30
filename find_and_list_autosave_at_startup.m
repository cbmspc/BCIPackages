% This script is meant to be called by matlabrc during startup.

tmp_saveslist = dir([feature('logdir') filesep 'matlabexitsave*.mat']);
if ~isempty(tmp_saveslist)
    tmp_lastsave = [feature('logdir') filesep tmp_saveslist(end).name];
else
    tmp_lastsave = '';
end
if exist(tmp_lastsave, 'file') && now - tmp_saveslist(end).datenum < 31 %#ok<TNOW1>
    clear tmp_result tmp_saveslist
    fprintf('\n\n\x25ca Restoring workspace from exitsave: ');
    try %#ok<TRYNC>
        load(tmp_lastsave); fprintf('Variables\x2713. ');
        if exist('tmp_cwd','var')
            cd(tmp_cwd); fprintf('Folder\x2713. ');
        end
        if exist('tmp_custompath','var')
            path(tmp_custompath); fprintf('Paths\x2713. ');
        end
    end
    try
        movefile(tmp_lastsave, getappdata(0,'g_matlabbaseworkspaceautosavefile'), 'f');
    catch
        try delete(tmp_lastsave); end %#ok<TRYNC>
    end
    clear tmp_*
    fprintf('Done. \n\x25c8 Please note your current directory: %s\n\n', cd);
else
    clear tmp_*
    % If there is no exit save, see if there was an autosave
    tmp_autosavelist = dir([feature('logdir') filesep 'matlabbaseautosave*.mat']);
    if length(tmp_autosavelist) > 1
        %fprintf('\x25c8 %i autosave files from previous MATLAB instances are available. Click <a href="matlab: uiopen(''%s'')">here</a> to choose and load into workspace.\n', length(tmp_autosavelist), [feature('logdir') filesep 'matlabbaseautosave*.mat']);
        fprintf('\x25c8 Found one or more autosave files, which may contain workspace variables lost due to a crash or reboot:\n');
        [~,tmp_i] = sort([tmp_autosavelist.datenum],'descend');
        tmp_autosavelist = tmp_autosavelist(tmp_i);
        tmp_k = 0;
        for tmp_i = 1:length(tmp_autosavelist)
            tmp_regtokens = regexp(tmp_autosavelist(tmp_i).name, '^matlabbaseautosave(\d+)\.mat$', 'tokens', 'once');
            if ~isempty(tmp_regtokens) && (tmp_autosavelist(tmp_i).bytes > 128 || tmp_autosavelist(tmp_i).bytes == 0)
                tmp_dt = datetime(str2double(tmp_regtokens{1}) / 1000, 'ConvertFrom', 'epochtime', 'TimeZone', 'Etc/UTC');
                tmp_dt.TimeZone = 'local';
                tmp_dt.Format = "eee dd-MMM-uuuu HH:mm:ss.SSS";
                tmp_bytes = sum([tmp_autosavelist(startsWith({tmp_autosavelist.name},['matlabbaseautosave' tmp_regtokens{1}])).bytes]);
                tmp_k = tmp_k + 1;
                fprintf('   %i.\t<a href="matlab: loadparts(''%s'');">%s</a> with %s bytes of workspace data. Click the timestamp to load.\n', tmp_k, [tmp_autosavelist(tmp_i).folder filesep tmp_autosavelist(tmp_i).name], tmp_dt, addThousandsCommaSeparators(tmp_bytes));
            end
        end
    end
    clear tmp_autosavelist tmp_i tmp_regtokens tmp_dt tmp_bytes tmp_k
end
fprintf('\n\n');
