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
        fprintf('You have saved work for %s, last saved at %s. Type "cd %s; loadwork;" to load and resume.\n', a(i).name, g(1).date, a(i).name);
    end
end


