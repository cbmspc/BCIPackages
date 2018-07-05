function hash = md5sum (DataOrFile)
if exist(DataOrFile, 'file')
    hash = md5sumfile(DataOrFile);
else
    Opt.Method = 'MD5';
    Opt.Format = 'hex';
    hash = datahash(DataOrFile, Opt);
end