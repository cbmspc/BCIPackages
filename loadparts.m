%
% filepath must be absolute full path to file ending in .mat
% Do not append the suffix such as --small
%

function loadparts(filepath)
verbose = true;
workspace = 'base';
if ~endsWith(filepath, '.mat')
    error('Filepath must end with .mat');
end
if ~isfile(filepath)
    error('File does not exist: %s', filepath);
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
dlist = [dlist; dir([pathname filesep filename fileext])];
[~,i] = sort([dlist.bytes],'descend');
dlist = dlist(i);
fprintf('Loading workspace variables from %s: ', filepath);
if verbose
    wb = waitbar(0, 'Loading variables', 'Name', ['loadparts: ' filename]);
    set(findobj(wb,'Interpreter','tex'),'Interpreter','none');
else
    wb = -1;
end
totalbytes = sum([dlist(i).bytes]);
bytessofar = 0;
for i = 1:length(dlist)
    if dlist(i).bytes > 0
        if ishandle(wb)
            desc = strrep(regexprep(dlist(i).name, '^.*--(\w+)-.*$', '$1'),'_','-');
            if isempty(desc) || endsWith(desc,'.mat')
                desc = 'smaller variables';
            end
            waitbar(bytessofar/totalbytes,wb,sprintf('Loading %s', desc));
            drawnow
        end
        evalin(workspace, sprintf('load(''%s'');',get_cached_filepath([dlist(i).folder filesep dlist(i).name])));
        bytessofar = bytessofar + dlist(i).bytes;
    end
end
if ishandle(wb)
    delete(wb);
end
fprintf('done.\n');
