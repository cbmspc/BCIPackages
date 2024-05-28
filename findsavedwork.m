function findsavedwork()
pause(0.1);
b = pwd;
w = getworkdir();
a = dir(b);
for i = 1:length(a)
    d = [b filesep a(i).name];
    f = [w filesep 'matlab-' crc32(lower(getusername)) '-' crc32(lower(d)) '-index.mat'];
    g = dir(f);
    if ~isempty(g)
        fprintf('You have saved work for %s, last saved at %s. <a href="matlab: clear; cd(''%s%s%s''); loadwork;">Click here to clear workspace and load this work</a>.\n', a(i).name, g(1).date, getbciprogramsdir(), filesep, a(i).name);
        fprintf('Note that you can have multiple types of saved data:\n  1. "Saved Work" (this feature) where you manually saved your workspace specific to the folder you are in\n  2. Exit Save (when you close MATLAB properly, it will reload the workspace when you start MATLAB next time)\n  3. Auto Save (every 15-30 minutes in the background, in case MATLAB crashes)\n');
    end
end


