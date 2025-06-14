function b = getfreespacebytes(Folder)
if ~exist('Folder', 'var') || (~ischar(Folder) && ~isstring(Folder)) || ~isfolder(Folder)
    Folder = feature('logdir');
end
b = java.io.File(Folder).getUsableSpace;

