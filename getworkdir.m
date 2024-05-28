% This is currently hardcoded to D:\Users\YourUserName\work
function D = getworkdir ()
D = gettmpdir();

% First, if the user has a Documents folder, return the MATLAB folder
% inside it as the working directory
UPF = getenv('USERPROFILE');
if exist([UPF filesep 'Documents'], 'dir')
    if ~exist([UPF filesep 'Documents' filesep 'MATLAB'], 'dir')
        try
            mkdir([UPF filesep 'Documents' filesep 'MATLAB']);
            D = [UPF filesep 'Documents' filesep 'MATLAB'];
            return
        catch
            D = gettmpdir();
            return
        end
    end
    D = [UPF filesep 'Documents' filesep 'MATLAB'];
    return
end

% Then, if the user has a OneDrive and it has a Documents folder, return
% the MATLAB folder inside it
ODF = getenv('ONEDRIVE');
if exist([ODF filesep 'Documents'], 'dir')
    if ~exist([ODF filesep 'Documents' filesep 'MATLAB'], 'dir')
        try
            mkdir([ODF filesep 'Documents' filesep 'MATLAB']);
            D = [ODF filesep 'Documents' filesep 'MATLAB'];
            return
        catch
            D = gettmpdir();
            return
        end
    end
    D = [ODF filesep 'Documents' filesep 'MATLAB'];
    return
end

