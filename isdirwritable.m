function b = isdirwritable (d)
b = 0;
if ~exist(d,'dir')
    b = 0;
    return
end
%Test mkdir
fp = [d filesep rand_az(1,8) '.tmp'];
success = mkdir(fp);
if ~success || ~exist(fp,'dir')
    b = 0;
    try
        rmdir(fp);
    catch
    end
    return
else
    %mkdir successful. Test rmdir
    try
        rmdir(fp);
    catch
    end
    if exist(fp,'dir')
        b = 0;
        return
    end
end

%Test create file
fp = [d filesep rand_az(1,8) '.tmp'];
try
    fid = fopen(fp,'w');
    fwrite(fid,0);
    fclose(fid);
catch
end
if ~exist(fp,'file')
    b = 0;
    try
        delete(fp);
    catch
    end
    return
else
    %Write successful. Test deletion
    try
        delete(fp);
    catch
    end
    if exist(fp,'file')
        b = 0;
        return
    end
end

%All tests passed
b = 1;
return