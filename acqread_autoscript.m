global gAcqRead_AutoScript_Last_Folder_Path
if ~isempty(gAcqRead_AutoScript_Last_Folder_Path)
    [tmp_f, gAcqRead_AutoScript_Last_Folder_Path] = uigetfile([gAcqRead_AutoScript_Last_Folder_Path '*.acq']);
else
    [tmp_f, gAcqRead_AutoScript_Last_Folder_Path] = uigetfile('*.acq');
end
BiopacAcqFile = [gAcqRead_AutoScript_Last_Folder_Path tmp_f];
[biopacdata, biopacsamplerate, biopacchannames] = acqreadwrapper(BiopacAcqFile);
