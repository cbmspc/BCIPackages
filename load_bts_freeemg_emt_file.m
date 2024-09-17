% ASCII Importer to load BTS Bioengineering EMT file for 
% FreeEMG wireless EMG sensor recordings
%
% Input argument emt_file: Enter the full path to the file, 
% e.g. 'D:\Research\EMG\Subject_555\emgrecord_0555_00001.emt'
%
% Output data will be in microvolts
%
% This reads from an ASCII file and is slow. Please be patient.
%
%
function [emgdata, Fs, emgchannames] = load_bts_freeemg_emt_file (emt_file)
fh = fopen(emt_file, 'r');

NumChan = 0;
NumSamples = 0; %#ok<NASGU>
ChanNames = {};
Data = [];
k = 0;


while ~feof(fh)
    line = fgetl(fh);
    clear b
    try
        b = regexp(line, '^Measure unit:\s+(\w+)', 'tokens', 'once');
    catch
        % File is blank
        emgdata = [];
        Fs = 0;
        emgchannames = {};
        return
    end
    if ~isempty(b)
        MeasureUnit = b{1};
        continue
    end
    
    b = regexp(line, '^Tracks:\s+(\d+)', 'tokens', 'once');
    if ~isempty(b)
        NumChan = str2double(b{1});
        ChanNames = cell(1,NumChan);
        continue
    end
    
    b = regexp(line, '^Frequency:\s+(\d+)\s+\[(\w+)\]', 'tokens', 'once');
    if ~isempty(b)
        SampleRateValue = str2double(b{1}); %#ok<NASGU>
        SampleRateUnit = b{2}; %#ok<NASGU>
        continue
    end
    
    b = regexp(line, '^Frames:\s+(\d+)', 'tokens', 'once');
    if ~isempty(b)
        NumSamples = str2double(b{1});
        Time = nan(NumSamples, 1);
        Data = nan(NumSamples, NumChan);
        continue
    end
    
    b = regexp(line, ['^\s+Frame\s+Time\s+' repmat('([^\t]+)\t', 1, NumChan) '$'], 'tokens', 'once');
    if ~isempty(b)
        ChanNames = b;
        continue
    end
    
    b = regexp(line, ['^\s+(\d+)\s+(\d+\.\d+)\s+' repmat('(-?\d+\.\d+)\s+', 1, NumChan-1) '(-?\d+\.\d+)$'], 'tokens', 'once');
    if ~isempty(b)
        k = k + 1;
        Time(k) = str2double(b{2});
        Data(k,:) = str2double(b(3:end)); %#ok<AGROW>
        
        
        continue
    end
    
    
end

fclose(fh);


switch MeasureUnit
    case 'mV'
        emgdata = Data .* 1000;
    case 'uV'
        emgdata = Data;
end

Fs = round(1/median(diff(Time)));
emgchannames = regexprep(regexprep(regexprep(regexprep(ChanNames, '^ +', ''), ' +$', ''), ' +', ' '), 'EMG Signal (\d+)', 'emg$1');
