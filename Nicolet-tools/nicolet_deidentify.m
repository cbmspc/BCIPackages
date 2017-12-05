function nicolet_deidentify (File)
if ~exist(File, 'file')
    error('File not found: %s', File);
end
fprintf('Now reading %s\n', File);
a = dir(File);
FileSize = a.bytes;
[FirstName, LastName, FirstPos, RawFirstName, RawLastName, MRN, RawMRN, DOB, RawDOB] = nicolet_extract_patient_name(File);

fprintf('Extraction complete.\n\n');

seeksize = 1048576;

tmp = string_to_cell(DOB,'/');
DOBmonth = str2double(tmp{1});
DOBday = str2double(tmp{2});
DOByear = str2double(tmp{3});
DOBunix_time = (datenum(DOByear,DOBmonth,DOBday) - datenum(1970,1,1))*86400;
DOBdbl = (DOBunix_time + 2209161600)/24/3600;
DOBregexpattern = convert_hex_group_to_regex_pattern(convert_to_double_memory_format(DOBdbl));

if DOBday == 0 && FirstName(end) == 'X' && LastName(end) == 'X' && length(strrep(unique(MRN(2:end)),'-','')) == 1
    fprintf('Skipping. This file is alrady de-identified.\n\n');
    return
end


fh = fopen(File, 'r+');

tmp = regexp(File, '(SJ\d+)@(\d{8})_(\d{6})', 'tokens', 'once');
if ~isempty(tmp)
    tmp = cell_to_string(tmp, ' - ');
else
    tmp = '';
end

buf = '';
Pos1List = [];
Pos2List = [];
Pos3List = [];
Pos4List = [];
fprintf('Locating..\n');
wh = waitbar(0, tmp);
fseek(fh, 0, 'bof');
while ~feof(fh)
    prevbuf = buf;
    buf = fread(fh, seeksize, 'uint8=>char')';
    b1s = regexp([prevbuf buf], [RawFirstName 0 0 RawLastName], 'start', 'once');
    b2s = regexp([prevbuf buf], RawMRN, 'start', 'once');
    b3s = regexp([prevbuf buf], RawDOB, 'start', 'once');
    b4s = regexp([prevbuf buf], DOBregexpattern, 'start', 'once');
    FilePos = ftell(fh);
    if isempty(prevbuf)
        n = 1;
    else
        n = 2;
    end
    if ~isempty(b1s)
        Pos1 = FilePos - seeksize*n + b1s - 1;
        if Pos1 < 1, keyboard; end
        Pos1List = [Pos1List Pos1];
    end
    if ~isempty(b2s)
        Pos2 = FilePos - seeksize*n + b2s - 1;
        if Pos2 < 1, keyboard; end
        Pos2List = [Pos2List Pos2];
    end
    if ~isempty(b3s)
        Pos3 = FilePos - seeksize*n + b3s - 1;
        if Pos3 < 1, keyboard; end
        Pos3List = [Pos3List Pos3];
    end
    if ~isempty(b4s)
        Pos4 = FilePos - seeksize*n + b4s - 1;
        if Pos4 < 1, keyboard; end
        Pos4List = [Pos4List Pos4];
    end
    waitbar(FilePos/FileSize, wh);
end
delete(wh);

Pos1List = unique(Pos1List);
Pos2List = unique(Pos2List);
Pos3List = unique(Pos3List);
Pos4List = unique(Pos4List);

fprintf('%s %s can be found at %s\n', FirstName, LastName, num2str(Pos1List));
fprintf('%s can be found at %s\n', MRN, num2str(Pos2List));
fprintf('%s can be found at %s\n', DOB, num2str(Pos3List));
fprintf('Encoded DOB can be found at %s\n', num2str(Pos4List));

MatchString = [RawFirstName 0 0 RawLastName];
MaskedString = regexprep([RawFirstName 0 0 RawLastName], '(?<=[A-Za-z]\x00)[A-Za-z]', 'X');

MatchString2 = RawMRN;
MaskedString2 = regexprep(RawMRN, '(?<=(?:\d|-)\x00)\d', '8');

MatchString3 = RawDOB;
tmp = string_to_cell(RawDOB, '/');
tmp(1:2) = regexprep(tmp(1:2), '\d(?=\x00)', '0');
MaskedString3 = cell_to_string(tmp, '/');

pause(0.5);

for i = 1:length(Pos1List)
    fprintf('Writing over names at position %i\n', Pos1List(i));
    fseek(fh, Pos1List(i), 'bof');
    fwrite(fh, MaskedString);
end

for i = 1:length(Pos2List)
    fprintf('Writing over the MRN at position %i\n', Pos2List(i));
    fseek(fh, Pos2List(i), 'bof');
    fwrite(fh, MaskedString2);
end

for i = 1:length(Pos3List)
    fprintf('Writing over the DOB at position %i\n', Pos3List(i));
    fseek(fh, Pos3List(i), 'bof');
    fwrite(fh, MaskedString3);
end


MaskedDOBunix_time = (datenum(DOByear,1,1) - datenum(1970,1,1))*86400;
MaskedDOBdbl = (MaskedDOBunix_time + 2209161600)/24/3600;

for i = 1:length(Pos4List)
    fprintf('Writing over the encoded DOB at position %i\n', Pos4List(i));
    fseek(fh, Pos4List(i), 'bof');
    fwrite(fh, MaskedDOBdbl, 'double');
end


%\x4d\x98\x00\x00\x00\x00\x00\x00

fclose(fh);
fprintf('\n');
