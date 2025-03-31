% Work in progress.
% Not ready for production use
% Po, 2025-03-25

function find_and_load_saved_work()
pause(0.1);
S = settings();
hx = S.matlab.desktop.currentfolder.History.PersonalValue;
count = 0;
finfou = [0 0];
fprintf('Finding workspace files in the folders you used recently...\n');

dlist = repmat(dir,0);

tmp_autosavelist = dir([feature('logdir') filesep 'matlabbaseautosave*.mat']);
if ~isempty(tmp_autosavelist)
    dlist = [dlist; tmp_autosavelist];
end

for i = 1:length(hx)
    dlist = [dlist; findsavedworkinfolder(hx{i})];
end
if isempty(dlist)
    fprintf('Did not find any save files.\n');
end

{dlist.name}'

function dlist = findsavedworkinfolder(folderPath)
dlist = repmat(dir,0);
if ~isfolder(folderPath)
    return
end

hashed_dpath = datahash([folderPath filesep getusername()], struct('Input','array','Format','hex','Method','MD5'));
wildcards_to_search = {
    [folderPath filesep getusername() '-workspace*.mat']
    [get_nascache_dir() filesep 'saved-' hashed_dpath '-workspace*.mat']
    [folderPath filesep 'matlab.mat']
};
for i = 1:length(wildcards_to_search)
    dlist = [dlist; dir(wildcards_to_search{i})];
end


