function findsavedwork()
pause(0.1);
S = settings();
hx = S.matlab.desktop.currentfolder.History.PersonalValue;
count = 0;
finfou = [0 0];
fprintf('Finding workspace files in the folders you used recently...\n');
for i = 1:length(hx)
    filename = [getusername() '-workspace.mat'];
    if ~isfile([hx{i} filesep filename])
        filename = 'matlab.mat';
    end
    if exist([hx{i} filesep filename], 'file')
        count = count + 1;
        finfo = dir([hx{i} filesep filename]);
        if ~isequal(finfou, [finfo.bytes finfo.datenum])
            finfou = [finfo.bytes finfo.datenum];
            fprintf('%i. "<a href="matlab: cd(''%s''); loadwork;">%s</a>" found in %s:\n   %s bytes, last saved: %s (% 3i days ago)   <a href="matlab: cd(''%s''); loadwork;">LOAD</a>   <a href="matlab: uigetfile(''%s'');">BROWSE</a>   <a href="matlab: delete(''%s'');fprintf(''Deleted #%i\\n'');">DELETE</a>\n\n', ...
                count, hx{i}, filename, ...
                hx{i}, addSpacesBeforeString(addThousandsCommaSeparators(finfo.bytes),14), finfo.date, floor(now - finfo.datenum), ...
                hx{i}, [hx{i} filesep '*.mat'], [hx{i} filesep regexprep(filename, 'workspace\.mat$', 'workspace*.mat')], count  ); %#ok<TNOW1>
        end
    end
end
if count == 0
    fprintf('Did not find any save files.\n');
end

