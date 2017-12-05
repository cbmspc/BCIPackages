function [FirstName, LastName, FirstPos, RawFirstName, RawLastName, MRN, RawMRN, DOB, RawDOB] = nicolet_extract_patient_name (File)
seeksize = 1048576;
fh = fopen(File, 'r');
buf = '';
while ~feof(fh)
    prevbuf = buf;
    buf = fread(fh, seeksize, 'uint8=>char')';
    b1 = regexp([prevbuf buf], '\x00\x00([A-Z ]\x00){2,}\x00\x00([A-Z ]\x00){2,}\x00\x00', 'tokens', 'once');
    bs = regexp([prevbuf buf], '\x00\x00([A-Z ]\x00){2,}\x00\x00([A-Z ]\x00){2,}\x00\x00', 'start', 'once');
    if isempty(prevbuf)
        n = 1;
    else
        n = 2;
    end
    if ~isempty(b1)
        RawFirstName = b1{1};
        RawLastName = b1{2};
        b1 = regexprep(b1, '([A-Z])\x00', '$1');
        FirstName = b1{1};
        LastName = b1{2};
        fprintf(' FN = %s\n', FirstName);
        fprintf(' LN = %s\n', LastName);
        FirstPos = ftell(fh) - seeksize*n + bs + 2; % Remember Matlab starts an index at 1
    end
    
    b2 = regexp([prevbuf buf], '(\d\x00){9}\x00\x00', 'tokens', 'once');
    if isempty(b2)
        b2 = regexp([prevbuf buf], '([0-9\-]\x00){9}\x00\x00', 'tokens', 'once');
    end
    if ~isempty(b2)
        RawMRN = b2{1};
        MRN = regexprep(RawMRN, '\x00', '');
        fprintf('MRN = %s\n', MRN);
    end
    
    b3 = regexp([prevbuf buf], '(\d\x00)+\/\x00(\d\x00)+\/\x00(\d\x00)+\x00\x00', 'match', 'once');
    if ~isempty(b3)
        RawDOB = b3(1:end-2);
        DOB = regexprep(RawDOB, '\x00', '');
        fprintf('DOB = %s\n', DOB);
    end
    
    if ~isempty(b1) && ~isempty(b2) && ~isempty(b3)
        break;
    end
    
end
fclose(fh);
