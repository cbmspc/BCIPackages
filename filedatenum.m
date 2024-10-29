function filedatenum = filedatenum(filepath)
list = dir(filepath);
if isempty(list)
    filedatenum = NaN;
    return
end
filedatenum = list(1).datenum;
