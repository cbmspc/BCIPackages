
BiopacAcqFile = 'C:/data/abc.acq';

% Load Biopac data
disp('Biopac: loading ..');
[info, data] = acqread(BiopacAcqFile);

disp('Biopac: sample rate ..');
biopacsamplerate = 1000 / info.dSampleTime;

disp('Biopac: chan names ..');
biopacchannames = info.szCommentText;

disp('Biopac: data ..');
biopacdata = single(zeros(info.lBufLength(1),info.nChannels));
for ch = 1:info.nChannels
    biopacdata(:,ch) = single(data{ch}) .* info.dAmplScale(ch);
end
clear info data ch
disp('Biopac: done');

disp('Converting to double precision ..');
biopacdata = double(biopacdata);

