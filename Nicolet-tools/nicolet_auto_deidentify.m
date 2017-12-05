% Enter a folder and it'll deidentify all .e files
function nicolet_auto_deidentify (Folder, YourIndex, TotalIndex)
a = dir([Folder filesep]);
for i = YourIndex:TotalIndex:length(a)
    if a(i).isdir && a(i).name(1) ~= '.'
        nicolet_auto_deidentify([Folder filesep a(i).name], 1, 1);
    elseif ~a(i).isdir && ~isempty(regexp(a(i).name, '\.e$', 'match', 'once'))
        try
            nicolet_deidentify([Folder filesep a(i).name]);
        catch me
            fprintf('Error: %s\nIn: %s\n', me.message, me.stack(1).file);
        end
    end
end
