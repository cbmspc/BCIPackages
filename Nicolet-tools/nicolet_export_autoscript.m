function nicolet_export_autoscript()
[FN, PN] = uigetfile('D:\Tracie\*.e', 'Locate the Nicolet .e file');
b = regexp(PN, '\\(SJ\d+@\d+_\d+)\\', 'tokens', 'once');
if ~isempty(b)
    Outfolder = ['D:\ResearchWorkspace\neurourology' filesep b{1}];
    if ~exist(Outfolder, 'dir')
        mkdir(Outfolder);
        nicolet_export([PN FN], Outfolder);
        nicolet_merge_channels(Outfolder);
    else
        fprintf('The output folder "%s" already exists. You may have already exported it. If not, delete it before trying to export.\n', Outfolder);
    end
end
