function ClaimSuccessful = claim_task_to_process(SaveDir, Task_ID_String, timeout_minutes, DoneFileMask)
ClaimSuccessful = false;
if ~exist('timeout_minutes','var') || numel(timeout_minutes) ~= 1 || ~isnumeric(timeout_minutes)
    timeout_minutes = 60;
end
if ~exist(SaveDir, 'dir')
    return
end

if ~isempty(dir([SaveDir filesep DoneFileMask]))
    return
end

my_instance_id = sanitizefilename(sanitizename(getcomputername()));

filepathmask = [SaveDir filesep '.processing_' sanitizefilename(sanitizename(Task_ID_String)) '_*.tmp'];
my_filename = ['.processing_' sanitizefilename(sanitizename(Task_ID_String)) '_' my_instance_id '.tmp'];
my_filepath = [SaveDir filesep '.processing_' sanitizefilename(sanitizename(Task_ID_String)) '_' my_instance_id '.tmp'];

if isempty(Task_ID_String) && timeout_minutes == 0
    % Special case: Wipe all claims
    filepathmask = [SaveDir filesep '.$inprog_*_*.tmp'];
    list = dir(filepathmask);
    for i = 1:length(list)
       try delete([SaveDir filesep list(i).name]); end %#ok<TRYNC>
    end
    return
end

% Otherwise, Task_ID_String cannot be empty
if isempty(Task_ID_String)
    return
end

list = dir(filepathmask);
list = list(~strcmpi({list.name}, my_filename));
if ~isempty(list) && now - max([list.datenum]) <= timeout_minutes/1440 %#ok<TNOW1>
    ClaimSuccessful = false; % Another worker is already processing this
    return
end

% Delete old claim files
remove_these = {list([list.datenum] > timeout_minutes/1440).name};
for i = 1:length(remove_these)
    try delete([SaveDir filesep remove_these{i}]); end %#ok<TRYNC>
end

% Now we claim this.
try %#ok<TRYNC>
    fh = fopen(my_filepath, 'w');
    fclose(fh); 
    ClaimSuccessful = true;
end
