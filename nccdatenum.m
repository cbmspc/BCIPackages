function [y, matlabdatenum] = nccdatenum (x)
if length(x) > 1 && x(1) == 'D'
    x = x(2:end);
end
    
if ischar(x) && length(x) >= 4 && length(x) <= 6
    y = nccdatenum_2(x);
    matlabdatenum = datenum(y);
elseif ischar(x) && length(x) == 8
    x = [str2double(x(1:4)), str2double(x(5:6)), str2double(x(7:8))];
    y = nccdatenum_1(x);
    matlabdatenum = datenum(x);
elseif length(x) == 1
    x = str2double(string_to_cell(datestr(x, 'yyyy-mm-dd'), '-'));
    y = nccdatenum_1(x);
    matlabdatenum = datenum(x);
else
    y = nccdatenum_1(x);
    matlabdatenum = datenum(x);
end


function d = nccdatenum_1 (yyyymmdd)
yyyy = yyyymmdd(1);
mm = yyyymmdd(2);
dd = yyyymmdd(3);
yyyy = mod(yyyy,100);
d = [num2str(yyyy) dec2base(mm * 100 + dd, 36)];

function yyyymmdd = nccdatenum_2 (d)
mmdd = base2dec(d(3:4), 36);
yyyy = str2double(['20' d(1:2)]);
dd = mod(mmdd,100);
mm = (mmdd-dd)/100;
yyyymmdd = [yyyy, mm, dd];
