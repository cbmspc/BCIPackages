% Main script for acquiring training EEG


% Master channels list: get_eeg_sensor_montage.m
%
%

close all
setmatlabtitle('Nexus EEG');
fprintf('\n\n\n');
tempdir = gettmpdir();
if exist([tempdir '\Parameters-autosave.mat'],'file')
    a = dir([tempdir '\Parameters-autosave.mat']);
    if (now - a.datenum)*86400 <= 600
        load([tempdir '\Parameters-autosave.mat']);
        fprintf('Loaded parameters from auto-recovery.\n');
        delete([tempdir '\Parameters-autosave.mat']);
        fprintf('\n\n');
        cd(cwd);
        clear a cwd;
    end
end

if ~exist('naa_tstart','var')
    naa_tstart = tic;
end

fprintf('This computer is: %s \n', getcomputername());

if ~exist('naa_in1','var') || toc(naa_tstart) > 28800
    [naa_in1, ~, Montage, NexusNumber] = nexus_prompt_select_montage ('nexus_acq_autoscript');
else
    fprintf('Current Nexus config is: %s\n', naa_in1);
end

if ~exist('naa_in2','var') || toc(naa_tstart) > 1800
    naa_in2 = getsubjectid();
else
    fprintf('Current Subject is: %s\n', naa_in2);
end

naa_tstart = tic;

if ~exist('NexusAcqSampleRate','var') || isempty(NexusAcqSampleRate)
    q = listdlg('ListString', {'128 Hz', '256 Hz', '512 Hz', '1024 Hz', '2048 Hz'}, 'SelectionMode', 'single', 'PromptString', 'Nexus sampling frequency', 'InitialValue', 2);
    switch q
        case 1
            NexusAcqSampleRate = 128;
        case 2
            NexusAcqSampleRate = 256;
        case 3
            NexusAcqSampleRate = 512;
        case 4
            NexusAcqSampleRate = 1024;
        case 5
            NexusAcqSampleRate = 2048;
    end
end

if ~exist('SWcar','var') || isempty(SWcar)
    SWcar = 1;
end

if ~exist('SWplot','var') || isempty(SWplot)
    SWplot = 1;
end

if ~exist('SWfilter','var') || isempty(SWfilter) || toc(naa_tstart) > 28800
    SWfilter = 1;
end
if ~exist('Fcutoff','var') || isempty(Fcutoff) || toc(naa_tstart) > 28800
    Fcutoff = NexusAcqSampleRate / 6.4;
end
if ~exist('Fcutoffhigh','var') || isempty(Fcutoffhigh) || toc(naa_tstart) > 28800
    Fcutoffhigh = 1;
end
if ~exist('SWautoscale','var') || isempty(SWautoscale) || toc(naa_tstart) > 28800
    SWautoscale = 0;
end
if ~exist('globalscale','var') || isempty(globalscale) || toc(naa_tstart) > 28800
    globalscale = 0.25;
end
clear SubjectName

if strcmpi(naa_in1, 'userdefined')
    Montage.userdefined = tmp_NexusChanNames.ChanNames;
end

ChanNames = Montage.(naa_in1);
ConfigName = naa_in1;
NexusNum = num2str(NexusNumber, '%.02i');


% switch naa_in1
%     case 1
%         ChanNames = nexus1;
%         ConfigName = 'Single Nexus (#1), Electrodes arranged as labeled';
%     case 2
%         ChanNames = nexus2;
%         ConfigName = 'Single Nexus (#2), Electrodes arranged as labeled';
%     case 3
%         ChanNames = synfi;
%         ConfigName = 'SynFi (#1 + #2), Electrodes arranged as labeled';
%     case 5
%         ChanNames = zoranoc2;
%         ConfigName = 'Single Nexus, Electrodes arranged as Zoran Optimistic Configuration #2';
%     case 6
%         ChanNames = string_to_cell(num2str(1:33),' ');
%         ConfigName = 'Single Nexus (#1), Electrodes arranged numerically';
%     case 7
%         ChanNames = string_to_cell([num2str(1:33,'A%i ') ' ' num2str(1:33,'B%i ')],' ');
%         ConfigName = 'SynFi (#1 + $2), Electrodes arranged numerically';
%     case 8
%         ChanNames = CustomChanNames;
%         ConfigName = 'SynFi (#1 + $2), Electrodes arranged manually';
%     case 9
%         ChanNames = acticap1;
%         ConfigName = 'ActiCap #1';
%     otherwise
%         clear naa_in1
%         error('Montage set error.');
% end
% 
% NexusNum = num2str(naa_in1,'%.2i');
% if naa_in1 == 3
%     NexusNum = '';
% end

if isempty(naa_in2)
    clear naa_in2;
    error('Subject Name is empty');
end

SubjectName = getsubjectid(naa_in2);
if isempty(SubjectName)
    error('SubjectName is empty');
end
naa_in2 = SubjectName;

aa = dir([getdesktopdir() filesep SubjectName '-*nexus*.mat']);
NumRecordedRuns = length(aa);
LastRecordedRun = 0;
for i = 1:NumRecordedRuns
    tmp = regexp(aa(i).name, '-(\d+)nexus','tokens','once');
    if ~isempty(tmp) && str2double(tmp{1}) > LastRecordedRun
        LastRecordedRun = str2double(tmp{1});
    end
end


fprintf('There are already %i recorded sessions. The last session number is %i\n', NumRecordedRuns, LastRecordedRun);


retry = 1;
while retry
    %20131001
    %naa_in3 = input('To record, enter session number, or enter 0 to check signals: ');
    naa_in3 = 0;
    RunNum = num2str(naa_in3,'%.2i');
    if isempty(naa_in3) || ~isnumeric(naa_in3) || round(naa_in3) ~= naa_in3 || naa_in3 < 0 || isinf(naa_in3) || naa_in3 > 99
        fprintf('Run Number must be integer between 1 and 99 (inclusive)\n');
        retry = 1;
        continue
    elseif exist([getdesktopdir() filesep SubjectName '-' RunNum 'nexus' NexusNum '.mat'], 'file')
        fprintf('Run already recorded. The last recorded session was %02i\n', LastRecordedRun);
        retry = 1;
        continue
    elseif naa_in3 ~= 0 && naa_in3 ~= LastRecordedRun+1
        fprintf('Warning: You seem to have skipped a session number. \n');
        fprintf('         This should be session %02i but you entered %02i\n', LastRecordedRun+1, naa_in3);
        uin = input('Enter Y to confirm, N to enter another session number: ','s');
        if strcmpi(uin,'y')
            retry = 0;
            break
        end
    else
        retry = 0;
        break
    end
end

SJprivacy = 1;
if ~verify_consent(SubjectName, 1)
    clear naa_in2
    error('You did not verify who the subject is. Run nexus_acq_autoscript again and make sure you verify the subject''s identity.');
end

if naa_in3 == 0
    fprintf('\n\n\n');
    disp('Signal check');
    %nexus_timeout_check
    %nexus_autoretry_check
    UpdateStatus('EEG signal check started.');
    statuscode = nexus_scope(ChanNames, NexusAcqSampleRate, 'SWplot', SWplot, 'SWfilter', SWfilter, 'SWcar', SWcar, 'Fcutoff', Fcutoff, 'Fcutoffhigh', Fcutoffhigh, 'SWautoscale', SWautoscale, 'SaveDir', getdesktopdir(), 'SubjectName', SubjectName, 'globalscale', globalscale);
    naa_tnexusactive = tic;
    if ~isempty(statuscode) && isnan(statuscode(1))
        nexus_acq_autoscript_restart
    end
    clear statuscode
    UpdateStatus('EEG signal check completed.');
    return
end

disp(['Subject: ' SubjectName]);
if ~isempty(NexusNum)
    disp(['  Nexus: ' NexusNum]);
else
    disp(['  Nexus: ' 'synfi']);
end
disp(['Run: ' RunNum]);

close all
close all hidden
nexus_timeout_check
fighand = figure(floor(rand*10000));
set(fighand, 'DockControls', 'off', 'MenuBar', 'none', 'NumberTitle', 'off');
set(gca, 'Xlim', [0 1], 'Ylim', [0 1]);
axis off
NexusNumDisp = NexusNum;
if isempty(NexusNum)
    NexusNumDisp = 'synfi';
end

%nexus_autoretry_check
clear rawdata;

if length(ChanNames) > 34
    % If using other than Nexus1, check memory 
    % If using Nexus1, possibly using Bluetooth, thus don't check.
    [MaxRecordTime, tmp] = nexus_scope(ChanNames, NexusAcqSampleRate, 'SWmemcheck', 1);
    MRTIME = [num2str(floor(MaxRecordTime/60)) ' minutes'];
else
    MRTIME = 'not checked';
    MaxRecordTime = -1;
end


thand(1) = text(0.5, 1.00, 'READY TO START', 'FontSize', 24, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'BackgroundColor', [0.9 0.9 0.9], 'Color', [0.1 0.5 0.1], 'FontWeight', 'bold');
thand(2) = text(0.0, 0.9, ['Subject: ' SubjectName 10 '  Nexus: ' NexusNumDisp ' (' num2str(length(ChanNames)) ' chans, ' num2str(NexusAcqSampleRate) ' Hz)' 10 'Run: ' RunNum 10 '   Date: ' datestr(now,'yyyy-mm-dd HH:MM')], 'FontSize', 16, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontName', 'FixedWidth');
thand(3) = text(0.10, 0.0, ['Press ENTER to start recording for real.' 10 'Or, press N and ENTER to cancel, R to restart.'], 'FontSize', 14, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'BackgroundColor', [0.95 0.95 0.95], 'Color', [0.1 0.4 0.1]);
thand(4) = text(0.5, 0.5, '', 'FontSize', 18, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'BackgroundColor', [0.9 0.9 0.1], 'Color', [0.9 0.1 0.1], 'FontWeight', 'bold');
thand(5) = text(0.5, 0.45, ['Maximum recording time: ' MRTIME], 'FontSize', 14, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

disp(['Estimated maximum recording time: ' MRTIME]);
if MaxRecordTime > 0 && floor(MaxRecordTime) < 25*60
    set(thand(4), 'String', '               WARNING               ');
    set(thand(5), 'String', ['Only ' num2str(floor(MaxRecordTime/60)) ' minutes data can be stored.' 10 '']);
end

disp(' ');
UpdateStatus('EEG recording ready to start.');
ticprompt = tic;
uin = input('Enter Y to start recording for real, N to cancel, R to restart: ','s');
if toc(ticprompt) > 120 || (~isempty(uin) && strcmpi(uin(1),'N'))
    try %#ok<TRYNC>
        figure(1);
        set(thand(1), 'String', 'Exiting program');
        for i = 2:5
            set(thand(i), 'String', '');
        end
        pause(1.0);
    end
    close all
    UpdateStatus('EEG recording canceled.');
    disp('Program exited.');
    return
elseif ~isempty(uin) && strcmpi(uin(1),'R')
    try %#ok<TRYNC>
        figure(1);
        set(thand(1), 'String', 'Restarting MATLAB');
        for i = 2:5
            set(thand(i), 'String', '');
        end
        pause(1.0);
    end
    disp('Restarting MATLAB.');
    nexus_acq_autoscript_restart
end
close all
clear rawdata ADCfactor
UpdateStatus('EEG recording started.');
[rawdata, NexusAcqSampleRate, ADCfactor, ~, ~, suid] = nexus_scope(ChanNames, NexusAcqSampleRate, 'SWplot', SWplot, 'SWfilter', SWfilter, 'SWcar', SWcar, 'SWautoscale', SWautoscale, 'SaveDir', getdesktopdir());
naa_tnexusactive = tic;

if isnan(rawdata(1))
    fprintf('\n');
    fprintf('\n');
    fprintf('An error occurred during acquisition, this session will be discarded.\n');
    fprintf('\n');
    fprintf('\n');
    nexus_acq_autoscript_restart
    return
end

sh = waitbar(0, 'Saving file...');
set(sh, 'CloseRequestFcn', '');
LocalSaveFile = [getdesktopdir() filesep SubjectName '-' RunNum 'nexus' NexusNum '.mat'];
save(LocalSaveFile,'rawdata','NexusAcqSampleRate','ADCfactor','synfi','nexus1','nexus2','ChanNames','SWcar','SWfilter','SWplot','zoranoc1','zoranoc2');
if ishandle(sh), waitbar(1/3,sh); end
md5_local = md5sum(LocalSaveFile);
fprintf('Saved as %s (md5sum = %s)\n', LocalSaveFile, md5_local);
if ishandle(sh), waitbar(2/3,sh); end

RemoteSavePath = 'D:\mnt\Desktop';
RemoteSaveFile = [RemoteSavePath filesep SubjectName '-' RunNum 'nexus' NexusNum '.mat'];
LockFile = 'Delete this when finished copying.tmp';
if exist([RemoteSavePath filesep LockFile],'file')
    fileattrib([RemoteSavePath filesep LockFile],'+w');
end
fh = fopen([RemoteSavePath filesep LockFile],'w');
fclose(fh);

fprintf('Copying experiment file to database computer... ');
suc = copyfileverify(LocalSaveFile, RemoteSaveFile);

if suc
    fprintf('done.\n');
    delete([RemoteSavePath filesep LockFile]);
    %fprintf('Copied to the Desktop of database computer\n');
else
    fprintf('copy failed.\n');
    fprintf('Failed to copy to the Desktop of database computer. Please do it manually.\n');
    uiwait(warndlg('Copying of the experiment file to the database computer has failed. You have to do it manually.', 'Copy failed'));
end
if ishandle(sh), waitbar(1,sh); end

UpdateStatus('EEG recording completed.');
if ishandle(sh), delete(sh); end

