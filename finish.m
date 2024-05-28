% Sourced from c:\bin\matlab_packages\finish.m
fprintf('MATLAB is exiting...');
evalin('base', 'tmp_cwd = pwd;');
evalin('base', 'tmp_custompath = path;');
try
    evalin('base', 'save([feature(''logdir'') filesep ''matlabexitsave'' num2str(java.lang.System.currentTimeMillis) ''.mat''], ''-v7.3'', ''-nocompression'');');
    autosave delete
catch
end
