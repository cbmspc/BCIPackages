function hash = sha1sum (data)
Opt.Method = 'SHA-1';
Opt.Format = 'hex';
hash = datahash(data, Opt);
