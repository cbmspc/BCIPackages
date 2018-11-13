% ASCII Importer to load BTS Bioengineering EMT file
function [Data, Time, ChanNames, MeasureUnit] = load_bts_emt (emt_file)
fh = fopen(emt_file, 'r');

NumChan = 0;
NumSamples = 0;
ChanNames = {};
Data = [];
k = 0;
wh = waitbar(0, 'Please wait');

while ~feof(fh)
    line = fgetl(fh);
    b = regexp(line, '^Measure unit:\s+(\w+)', 'tokens', 'once');
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
        SampleRateValue = str2double(b{1});
        SampleRateUnit = b{2};
        continue
    end
    
    b = regexp(line, '^Frames:\s+(\d+)', 'tokens', 'once');
    if ~isempty(b)
        NumSamples = str2double(b{1});
        Time = nan(NumSamples, 1);
        Data = nan(NumSamples, NumChan);
        fprintf('Data resized to %i x %i.\n', NumSamples, NumChan);
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
        Data(k,:) = str2double(b(3:end));
        
        waitbar(k/NumSamples, wh);
        
        continue
    end
    
    
end

delete(wh);

fclose(fh);

