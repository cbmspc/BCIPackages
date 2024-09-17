function findsavedwork()
pause(0.1);
S = settings();
hx = S.matlab.desktop.currentfolder.History.PersonalValue;
count = 0;
finfou = [0 0];
for i = 1:length(hx)
    if exist(hx{i}, 'dir') && exist([hx{i} filesep 'matlab.mat'], 'file')
        count = count + 1;
        finfo = dir([hx{i} filesep 'matlab.mat']);
        if ~isequal(finfou, [finfo.bytes finfo.datenum])
            finfou = [finfo.bytes finfo.datenum];
            vinfo = whos('-file', [hx{i} filesep 'matlab.mat']);
            [~, ind] = sort([vinfo.bytes], 'descend');
            vbignames = {vinfo(ind).name};
            if length(vbignames) > 5
                vbignames = [cell_to_string(vbignames(1:5), ', ') ', ...'];
            else
                vbignames = cell_to_string(vbignames, ', ');
            end

            fprintf('▶ "%s" contains workspace variables you saved earlier.\n    Size: %i bytes\n    Count: %i variables\n    Last saved: %s\n    Variable names: %s\n    Click <a href="matlab: clear">here</a> to clear workspace, and then <a href="matlab: cd(''%s''); load;">here</a> to load this work.\n', hx{i}, sum([vinfo.bytes]), length(vinfo), finfo.date, vbignames, hx{i});
        end
    end
end

