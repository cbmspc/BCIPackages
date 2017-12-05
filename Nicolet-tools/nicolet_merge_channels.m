function nicolet_merge_channels (Folder)
persistent LastPath
if ~exist('Folder','var') || isempty(Folder) || ~exist(Folder,'dir')
    if ~isempty(LastPath) && exist(LastPath,'dir')
        PN = uigetdir(LastPath, 'Locate the base folder of the exported data');
    else
        PN = uigetdir('', 'Locate the base folder of the exported data');
    end
    if PN ~= 0
        Folder = PN;
        LastPath = PN;
    else
        error('User canceled.');
    end
end

a_test = dir([Folder filesep 'Segment_1']);
a = dir([Folder filesep 'Segment_*']);
if isempty(a_test)
    error('Cannot find any segments. Make sure you selected the base folder, and not a segment subfolder.');
end


Nsegment = length(a);
fprintf('There are %i segments.\n', Nsegment);

for seg = 1:Nsegment
    if ~a(seg).isdir
        continue
    end
    Segmentfolder = [Folder filesep a(seg).name];
    S = load([Segmentfolder filesep 'Segment_info.mat']);
    Nchan = length(S.ChannelNames);
    fprintf('%s started at %s, has %i channels, and is %g seconds long.\n', a(seg).name, S.StartDateStr, Nchan, S.DurationSeconds);
    
    SampleRates = zeros(1,Nchan);
    ChannelNames = cell(1,Nchan);
    Scales = zeros(1,Nchan);
    for ch = 1:Nchan
        if exist([Segmentfolder filesep sprintf('Channel_%i_info.mat', ch)], 'file')
            S = load([Segmentfolder filesep sprintf('Channel_%i_info.mat', ch)]);
            ChannelNames{ch} = S.ChannelName;
            SampleRates(ch) = S.SamplesPerSecond;
            Scales(ch) = S.Scale;
        end
    end
    UniqueSampleRates = unique(SampleRates);
    if length(UniqueSampleRates) == 1 && UniqueSampleRates == 0
        fprintf('There are no channels to merge in this segment.\n');
        continue
    end
    
    UniqueSampleRates = UniqueSampleRates(UniqueSampleRates>0);
    
    Nusr = length(UniqueSampleRates);
    fprintf('There are %i unique sample rates among the channels.\n', Nusr);
    
    for u = 1:Nusr
        
        SampleRate = UniqueSampleRates(u);
        ChanNames = ChannelNames(SampleRate==SampleRates); %#ok<NASGU>
        ChanScales = Scales(SampleRate==SampleRates); %#ok<NASGU>
        IndexList = find(SampleRate==SampleRates);
        fprintf('Now combining all channels with Sample Rate = %g ..\n', SampleRate);
        clear data
        ch = IndexList(1);
        S = load([Segmentfolder filesep sprintf('Channel_%i_data.mat', ch)]);
        Datalength = numel(S.data);
        data = zeros(Datalength, length(IndexList));
        data(:,1) = S.data;
        for k = 2:length(IndexList)
            ch = IndexList(k);
            S = load([Segmentfolder filesep sprintf('Channel_%i_data.mat', ch)]);
            data(:,k) = S.data;
        end
        fprintf('Saving into one file..\n');
        save([Segmentfolder filesep sprintf('Dataset_%i_data.mat', u)], 'data', '-v7.3');
        save([Segmentfolder filesep sprintf('Dataset_%i_info.mat', u)], 'SampleRate', 'ChanNames', 'ChanScales', '-v7.3');
        
        fprintf('Deleting redundant unmerged files..\n');
        for k = 1:length(IndexList)
            ch = IndexList(k);
            retry = 1;
            while retry
                try
                    delete([Segmentfolder filesep sprintf('Channel_%i_data.mat', ch)]);
                    delete([Segmentfolder filesep sprintf('Channel_%i_info.mat', ch)]);
                    retry = 0;
                catch
                    fprintf('Could not delete %s. Retrying..\n', [Segmentfolder filesep sprintf('Channel_%i_*.mat', ch)]);
                    pause(5.0);
                end
            end
        end
    end
    
    
    
    
end

