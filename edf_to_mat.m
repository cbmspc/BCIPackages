% Mainly use this for Natus edf exports
%
function edf_to_mat(edffilepath)
persistent lastpath
if nargin == 0
    [f, p] = uigetfile('*.edf', 'Select an EDF file', lastpath);
    edffilepath = [p f];
end
savefilepath = [p regexprep(f, '\.edf$', '') '.mat'];
if isfile(savefilepath)
    error('%s already exists.', savefilepath);
end
fprintf('Loading the EDF file %s\n', edffilepath);
[hdr, record] = edfread(edffilepath);
fprintf('Reformatting... ');
c = string_to_cell(hdr.patientID,' ');
l = true(1,length(c));
mrn = '';
sex = '';
age = 0;
for i = 1:length(c)
    b0 = regexp(c{i}, '^[0-9]+$', 'match', 'once');
    if ~isempty(b0) && l(i) && isempty(mrn)
        % Medical record number
        mrn = regexprep(b0, '[0-9]', '*');
        l(i) = false;
        continue
    end
    if l(i) && isempty(sex)
        b1 = regexpi(c{i}, '^M$|^F$|^Male|^Female', 'match', 'once');
        if ~isempty(b1)
            sex = b1(1);
            l(i) = false;
            continue
        else
            sex = 'U';
        end
    end
    if l(i) && age == 0
        try
            b2 = datetime(c{i});
            % Age
            age = floor(years(datetime("now") - b2));
            l(i) = false;
            continue
        catch
            age = 0;
        end
    end
end
name = cell_to_string(c(l), ' ');
if contains(name,',')
    name = regexprep(name, '^([^,]+),([^,]+)$', '$2 $1');
end
name = regexprep(regexprep(name, '^\s+', ''), '\s+$', '');
name = strrep(name, '-', ' ');
name = strrep(name, '_', ' ');
name = strrep(regexprep(name, '(\S)(\S+)', '$1'), ' ', '');
hdr.patientID = [sex ' ' num2str(age) ' ' name];
hdr.recordID = '';
SubjectInfo = hdr.patientID;
data = record.'; clear record
Fs = size(data,1) / (hdr.records * hdr.duration);
ChanNames = hdr.label;
fprintf('Saving to disk... ');
form = '-v7.3';
w = whos('data');
if w.bytes < 2^31
    form = '-v7';
end
save(savefilepath, form, 'data', 'Fs', 'ChanNames', 'SubjectInfo');
fprintf('done. Save file: %s\n', savefilepath);
lastpath = p;
