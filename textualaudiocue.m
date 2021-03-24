function textualaudiocue (Ntrial, IdleDuration, WalkDuration, SWheadphone, IdleString, WalkString)
% Gives the cues IDLE, WALK for a total number of Ntrial times.
% Also produces sound for re-acquisition.

setmatlabtitle('CUE');

ProgramStart = tic;

%Idle.String = '';
Idle.Color = [0.9 0.1 0.1];
Idle.FontSize = 96;
Idle.FontName = 'Times New Roman';
%Walk.String = '';
Walk.Color = [0.1 0.9 0.1];
Walk.FontSize = 96;
Walk.FontName = 'Times New Roman';

if nargin < 3
    error('Requires 3 arguments.');
end

if ~exist('SWheadphone','var') || isempty(SWheadphone)
    SWheadphone = 0;
end

if ~exist('IdleString','var') || isempty(IdleString)
    IdleString = 'IDLE';
end

if ~exist('WalkString','var') || isempty(WalkString)
    WalkString = 'WALK';
end

Idle.String = IdleString;
Walk.String = WalkString;

%% Set up audio signals
if SWheadphone
    AStruct = load('VoiceStartStop.mat');
    Audio.Fs = AStruct.VoiceFs;
    Audio.IdleY = AStruct.VoiceStop * 1e-1;
    Audio.WalkY = AStruct.VoiceStart * 1e-1;
else
    Audio.Lag = 0.03;
    Audio.Fs = 22050;
    Audio.IdleFreq = 220;
    Audio.IdleAmp = 0;
    Audio.WalkFreq = 350;
    Audio.WalkAmp = 0.750;
    Audio.IdleT = (1:Audio.Fs*(IdleDuration-Audio.Lag))/Audio.Fs;
    Audio.WalkT = (1:Audio.Fs*(WalkDuration-Audio.Lag))/Audio.Fs;
    Audio.IdleY = Audio.IdleAmp*sin(2*pi*Audio.IdleFreq*Audio.IdleT);
    Audio.WalkY = Audio.WalkAmp*sin(2*pi*Audio.WalkFreq*Audio.WalkT);
end

%% Set up the User Figure
close all
MonitorPositions = get(0,'MonitorPositions');
if size(MonitorPositions,1) < 2
    warning('No second monitor');
else
    UserFigurePosition = MonitorPositions(2,:) + [-1 1 1-MonitorPositions(1,3) 0];
    UserFigurePosition(2) = MonitorPositions(1,4) - MonitorPositions(2,4);
end

figure(1);
set(gcf,'Toolbar','none','MenuBar','none', 'NumberTitle','off','DoubleBuffer','on');
try
    set(gcf,'Position',UserFigurePosition);
end
hold on
set(gca,'Position',[0 0 1 1],'XTick',[],'YTick',[],'XLim',[-1 1],'YLim',[-1 1], 'Color',[0 0 0]);
sbuffer = text(0,0,'','HorizontalAlignment','center');
s2buffer = text(0.0,-0.2,'','HorizontalAlignment','center','FontSize',24,'Color',[0 0.5 0.5]);
fprintf('[%.3f] User Figure window configured. Get ready to start ..\n',toc(ProgramStart));
%UpdateInstantStatus('started', 'clearall');

q = questdlg('Ready to start texual cue?', 'textualaudiocue');
drawnow
switch q
    case 'Yes'
    otherwise
        close all
        return
end
    
EstTime = Ntrial * (IdleDuration + WalkDuration)/2 + 3.0;
%UpdateCueStatus(0, EstTime);   %JL 3/10/21
pause(3.0);

CuesGivenC = cell(1,Ntrial);

%% Start



for tri = 1:Ntrial
    cui = mod(tri-1,2)+1;
    switch cui
        case 1
            CueDuration = IdleDuration;
            CueString = Idle.String;
            Props = struct2cellpairs(Idle);
        case 2
            CueDuration = WalkDuration;
            CueString = Walk.String;
            Props = struct2cellpairs(Walk);
    end
    if length(CueDuration) > 1 && CueDuration(2) > CueDuration(1)
        CueDuration = CueDuration(1) + rand*(CueDuration(2) - CueDuration(1));
    else
        CueDuration = CueDuration(1);
    end
    CueGivenTime = toc(ProgramStart);
    fprintf('[%.3f] Trial #%3i, Cue#%3i, CueDuration=%6.3f, CueText: %s\n',CueGivenTime, tri, cui, CueDuration, CueString);
    UpdateInstantStatus(['TRIAL_' num2str(tri,'%03i')], ['TRIAL_' num2str(tri-1,'%03i')]);
    set(sbuffer, Props{:});
    Extent1 = get(sbuffer,'Extent');
    Extent2 = get(s2buffer,'Extent');
    bottom = Extent1(2);
    Pos2 = get(s2buffer,'Position');
    Pos2(2) = [bottom - Extent2(4)];
    set(s2buffer, 'String', [num2str(Ntrial-tri)], 'Position', Pos2);
    CuesGivenC{tri} = {CueGivenTime, CueDuration, CueString};
    drawnow
    audioplaystart = tic;
    switch cui
        case 1
            wavplay(Audio.IdleY,Audio.Fs);
        case 2
            wavplay(Audio.WalkY,Audio.Fs);
    end
    %CueDuration - toc(audioplaystart)
    pause(CueDuration - toc(audioplaystart));
end

Props = struct2cellpairs(Idle);
set(sbuffer, Props{:});
Extent1 = get(sbuffer,'Extent');
Extent2 = get(s2buffer,'Extent');
bottom = Extent1(2);
Pos2 = get(s2buffer,'Position');
Pos2(2) = [bottom - Extent2(4)];
set(s2buffer, 'String', 'End of session.', 'Position', Pos2);
pause(3.5);
fprintf('[%.3f] End of session.\n',toc(ProgramStart));
UpdateInstantStatus('done', 'clearall');
close all
