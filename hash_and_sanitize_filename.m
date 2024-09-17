function outputname = hash_and_sanitize_filename (inputname, truncatelength)
outputname = [crc32(lower(inputname)) '-' sanitizefilename(inputname)];
if ~exist('truncatelength', 'var') || isempty(truncatelength) || ~isnumeric(truncatelength) || ~isfinite(truncatelength)
    truncatelength = 64;
end
if truncatelength < 8
    truncatelength = 8;
end
if length(outputname) > truncatelength
    outputname = outputname(1:truncatelength);
end
