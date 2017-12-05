% This exports the Nicolet .e file into Matlab .mat files in subfolders
%
% written by Po T Wang
% uses Nicolet reader by jwagenaar at https://github.com/ieeg-portal/Nicolet-Reader

% OutFolder
% --- Segment 1
% ------ Segment metadata
% ------ Channel 1 data
% ------ Channel 2 data
% ------ Channel 3 data
% ------ Channel ... data
% --- Segment 2
% ------ Segment metadata
% ------ Channel 1 data
% ------ Channel 2 data
% ------ Channel 3 data
% ------ Channel ... data
% ...




function nicolet_export (File, OutFolder)
persistent LastPath
if ~exist('File','var') || isempty(File) || ~exist(File,'file')
    if ~isempty(LastPath) && exist(LastPath,'dir')
        [FN, PN] = uigetfile([LastPath '*.e'], 'Locate the Nicolet .e file');
    else
        [FN, PN] = uigetfile('*.e', 'Locate the Nicolet .e file');
    end
    if FN ~= 0
        File = [PN FN];
        LastPath = PN;
    else
        error('User canceled.');
    end
end

if ~exist('OutFolder','var') || isempty(OutFolder) || ~exist(OutFolder,'dir')
    if ~isempty(LastPath) && exist(LastPath,'dir')
        PN = uigetdir(LastPath, ['Locate a folder in which to save the exported data for ' File]);
    else
        PN = uigetdir('', ['Locate a folder in which to save the exported data for ' File]);
    end
    if PN ~= 0
        OutFolder = PN;
    else
        error('User canceled.');
    end
end


OBJ = NicoletFile(File);
Nsegment = length(OBJ.segments);
rootsavedir = OutFolder;
if ~exist(rootsavedir, 'dir')
    mkdir(rootsavedir);
end
for seg = 1:Nsegment
    % Note: Each segment can have completely different set of channels
    % Note: Each channel can have a different sample rate.
    savedir = [OutFolder filesep sprintf('Segment_%i', seg)];
    if ~exist(savedir, 'dir')
        mkdir(savedir);
    end
    SourceFileName = OBJ.fileName;
    StartDateNum = datenum(OBJ.segments(seg).dateStr);
    StartDateStr = OBJ.segments(seg).dateStr; %#ok<NASGU>
    DurationSeconds = OBJ.segments(seg).duration;
    ChannelNames = OBJ.segments(seg).chName;
    save([savedir filesep 'Segment_info.mat'], 'SourceFileName', 'StartDateNum', 'StartDateStr', 'DurationSeconds', 'ChannelNames', '-v7.3');
    
    for ch = 1:length(ChannelNames)
        ChannelName = OBJ.segments(seg).chName{ch};
        Scale = OBJ.segments(seg).scale(ch);
        SamplesPerSecond = OBJ.segments(seg).samplingRate(ch);
        
        if exist([savedir filesep sprintf('Channel_%i_info', ch) '.mat'], 'file') && exist([savedir filesep sprintf('Channel_%i_data', ch) '.mat'], 'file')
            % Files already exist. Check if all info are identical
            S = load([savedir filesep sprintf('Channel_%i_info', ch) '.mat']);
            Identical = 0;
            if isfield(S,'SourceFileName') && isfield(S,'StartDateNum') && isfield(S,'ChannelName') && isfield(S,'SamplesPerSecond') && isfield(S,'DurationSeconds') && isfield(S,'Scale')
                if ischar(S.SourceFileName) && strcmp(S.SourceFileName, SourceFileName) && ...
                        isnumeric(S.StartDateNum) && S.StartDateNum == StartDateNum && ...
                        ischar(S.ChannelName) && strcmp(S.ChannelName,ChannelName) && ...
                        isnumeric(S.SamplesPerSecond) && S.SamplesPerSecond == SamplesPerSecond && ...
                        isnumeric(S.DurationSeconds) && S.DurationSeconds == DurationSeconds && ...
                        isnumeric(S.Scale) && S.Scale == Scale
                    Identical = 1;
                end
            end
            if Identical
                fprintf('Skipping segment %i channel %i (already exported)\n', seg, ch);
                continue
            end
        end
        try
            data = getdata(OBJ, seg, [1 DurationSeconds*SamplesPerSecond], ch); %#ok<NASGU>
            % Save right away
            fprintf('Saving data for segment %i channel %i ..\n', seg, ch);
            save([savedir filesep sprintf('Channel_%i_data', ch) '.mat'], 'data', '-v7.3');
            save([savedir filesep sprintf('Channel_%i_info', ch) '.mat'], 'ChannelName', 'SamplesPerSecond', 'StartDateNum', 'StartDateStr', 'DurationSeconds', 'Scale', 'SourceFileName', '-v7.3');
        catch me
            fprintf('\nThere was an error when exporting %s, segment #%i, channel #%i (%s).\nThe program %s at line %i crashed in the function %s.\nError: %s (%s)\nThis affected channel segment was NOT exported.\n\n', File, seg, ch, ChannelName, me.stack(1).file, me.stack(1).line, me.stack(1).name, me.message, me.identifier);
        end
    end
end

return