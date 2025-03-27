%
% filepath must be absolute full path to file ending in .mat
% Do not append the suffix such as --small
%

function loadparts(filepath)
workspace = 'base';
if ~endsWith(filepath, '.mat')
    error('Filepath must end with .mat');
end
[pathname, filename, fileext] = fileparts(filepath);
if ~isfolder(pathname)
    error('Not a folder: %s', pathname);
end
if ~strcmp(fileext,'.mat')
    error('Filepath must end with .mat');
end

filename = regexprep(filename, '--small', '');

dlist = dir([pathname filesep filename '--*' fileext]);
[~,i] = sort([dlist.bytes],'descend');
dlist = dlist(i);
for i = 1:length(dlist)
    evalin(workspace, sprintf('load(''%s'');',get_cached_filepath([dlist(i).folder filesep dlist(i).name])));
end

