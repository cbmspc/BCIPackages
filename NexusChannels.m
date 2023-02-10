% Channel names GUI
%
%
% Master channels list: get_eeg_sensor_montage.m
%
%

function varargout = NexusChannels(varargin)
% NEXUSCHANNELS MATLAB code for NexusChannels.fig
%      NEXUSCHANNELS, by itself, creates a new NEXUSCHANNELS or raises the existing
%      singleton*.
%
%      H = NEXUSCHANNELS returns the handle to a new NEXUSCHANNELS or the handle to
%      the existing singleton*.
%
%      NEXUSCHANNELS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NEXUSCHANNELS.M with the given input arguments.
%
%      NEXUSCHANNELS('Property','Value',...) creates a new NEXUSCHANNELS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NexusChannels_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NexusChannels_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help NexusChannels

% Last Modified by GUIDE v2.5 20-Nov-2013 12:47:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NexusChannels_OpeningFcn, ...
                   'gui_OutputFcn',  @NexusChannels_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before NexusChannels is made visible.
function NexusChannels_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NexusChannels (see VARARGIN)

% Choose default command line output for NexusChannels
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes NexusChannels wait for user response (see UIRESUME)
% uiwait(handles.NexusChannels);

blankchannelmode_Callback(handles.blankchannelmode, 'init', handles);
ampsused_Callback(handles.ampsused, 'init', handles);

NexusChanNamesFile = [getdesktopdir() filesep 'NexusChanNames.mat'];
if exist(NexusChanNamesFile, 'file')
    tmp = load(NexusChanNamesFile, 'NexusChanNames');
    loadallchans(tmp.NexusChanNames, handles);
end



% --- Outputs from this function are returned to the command line.
function varargout = NexusChannels_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function ch101_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
updateall(hObject, eventdata, handles);

function ch102_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch103_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch104_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch105_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch106_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch107_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch108_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch109_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch110_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch111_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch112_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch113_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch114_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch115_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch116_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch117_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch118_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch119_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch120_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch121_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch122_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch123_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch124_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch125_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch126_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch127_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch128_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch129_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch130_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch131_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch132_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch201_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch202_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch203_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch204_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch205_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch206_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch207_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch208_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch209_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch210_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch211_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch212_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch213_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch214_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch215_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch216_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch217_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch218_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch219_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch220_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch221_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch222_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch223_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch224_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch225_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch226_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch227_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch228_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch229_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch230_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch231_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch232_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch1gnd_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch2gnd_Callback(hObject, eventdata, handles)
updateall(hObject, eventdata, handles);

function ch101_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>
editcfun(hObject, eventdata, handles);

function ch102_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch103_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch104_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch105_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch106_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch107_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch108_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch109_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch110_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch111_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch112_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch113_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch114_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch115_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch116_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch117_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch118_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch119_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch120_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch121_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch122_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch123_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch124_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch125_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch126_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch127_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch128_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch129_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch130_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch131_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch132_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch201_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch202_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch203_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch204_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch205_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch206_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch207_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch208_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch209_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch210_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch211_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch212_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch213_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch214_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch215_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch216_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch217_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch218_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch219_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch220_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch221_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch222_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch223_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch224_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch225_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch226_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch227_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch228_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch229_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch230_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch231_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch232_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch1gnd_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

function ch2gnd_CreateFcn(hObject, eventdata, handles)
editcfun(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1
imshow('nx32b.jpg', 'Parent', hObject);


% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2
imshow('nx32b.jpg', 'Parent', hObject);

function updateall(hObject, eventdata, handles)

amp = 1; ch = 0;
s = get(handles.ch1gnd, 'String');
s = regexp(s, '[A-Za-z0-9]+', 'match', 'once');
if isempty(s)
    bm = getappdata(handles.NexusChannels, 'BlankMode');
    if bm == 1
        s = ['NC' num2str(amp,'%i') '-' num2str(ch,'%i')];
    elseif bm == 2
        s = [char(64+amp) num2str(ch,'%i')];
    end
end
Gnd1 = s;

amp = 2; ch = 0;
s = get(handles.ch2gnd, 'String');
s = regexp(s, '[A-Za-z0-9]+', 'match', 'once');
if isempty(s)
    bm = getappdata(handles.NexusChannels, 'BlankMode');
    if bm == 1
        s = ['NC' num2str(amp,'%i') '-' num2str(ch,'%i')];
    elseif bm == 2
        s = [char(64+amp) num2str(ch,'%i')];
    end
end
Gnd2 = s;

for amp = 1:2
    channames = cell(1,32);
    for ch = 1:32
        h = handles.(['ch' num2str(amp,'%i') num2str(ch,'%02i')]);
        s = get(h, 'String');
        s = regexp(s, '[A-Za-z0-9]+', 'match', 'once');
        if isempty(s)
            bm = getappdata(handles.NexusChannels, 'BlankMode');
            if bm == 1
                s = ['NC' num2str(amp,'%i') '-' num2str(ch,'%i')];
            elseif bm == 2
                s = [char(64+amp) num2str(ch,'%i')];
            end
        end
        channames{ch} = s;
    end
    NexusChanNames.(['nexus' num2str(amp,'%i')]) = channames;
    if amp == 1
        NexusChanNames.synfi(1:32) = channames;
    elseif amp == 2
        NexusChanNames.synfi(34:65) = channames;
    end
end
NexusChanNames.nexus1{33} = 'D1';
NexusChanNames.nexus2{33} = 'D2';
NexusChanNames.nexus1{34} = 'DIAG';
NexusChanNames.nexus2{34} = 'DIAG';
NexusChanNames.synfi{33} = 'D1';
NexusChanNames.synfi{66} = 'D2';
NexusChanNames.synfi{67} = 'DIAG';
NexusChanNames.Gnd1 = Gnd1;
NexusChanNames.Gnd2 = Gnd2;
aum = getappdata(handles.NexusChannels, 'AmpsUsedMode');
switch aum
    case 1
        NexusChanNames.ChanNames = NexusChanNames.synfi;
    case 2
        NexusChanNames.ChanNames = NexusChanNames.nexus1;
    case 3
        NexusChanNames.ChanNames = NexusChanNames.nexus2;
end
assignin('base', 'tmp_NexusChanNames', NexusChanNames);
save([getdesktopdir() filesep 'NexusChanNames.mat'], 'NexusChanNames');


function editcfun(hObject, eventdata, handles)
set(hObject, 'String', '');
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function loadallchans (NexusChanNames, handles)
for ch = 1:32
    h = handles.(['ch1' num2str(ch,'%02i')]);
    s = NexusChanNames.nexus1{ch};
    e1 = ['NC1-' num2str(ch,'%i')];
    e2 = ['A' num2str(ch,'%i')];
    if strcmp(s, e1) || strcmp(s, e2)
        s = '';
    end
    set(h, 'String', s);
end
for ch = 1:32
    h = handles.(['ch2' num2str(ch,'%02i')]);
    s = NexusChanNames.nexus2{ch};
    e1 = ['NC2-' num2str(ch,'%i')];
    e2 = ['B' num2str(ch,'%i')];
    if strcmp(s, e1) || strcmp(s, e2)
        s = '';
    end
    set(h, 'String', s);
end
if isfield(NexusChanNames, 'Gnd1')
    set(handles.ch1gnd, 'String', NexusChanNames.Gnd1);
end
if isfield(NexusChanNames, 'Gnd2')
    set(handles.ch2gnd, 'String', NexusChanNames.Gnd2);
end

% --- Executes on selection change in ampsused.
function ampsused_Callback(hObject, eventdata, handles)
% hObject    handle to ampsused (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ampsused contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ampsused
sel = get(hObject,'Value');
switch sel
    case 1
        setappdata(handles.NexusChannels, 'AmpsUsedMode', 1);
    case 2
        setappdata(handles.NexusChannels, 'AmpsUsedMode', 2);
    case 3
        setappdata(handles.NexusChannels, 'AmpsUsedMode', 3);
end

if strcmp(eventdata, 'init')
    return;
end

updateall(hObject, eventdata, handles);



% --- Executes during object creation, after setting all properties.
function ampsused_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ampsused (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in blankchannelmode.
function blankchannelmode_Callback(hObject, eventdata, handles)
% hObject    handle to blankchannelmode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns blankchannelmode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from blankchannelmode
sel = get(hObject,'Value');
switch sel
    case 1
        setappdata(handles.NexusChannels, 'BlankMode', 1);
    case 2
        setappdata(handles.NexusChannels, 'BlankMode', 2);
end

if strcmp(eventdata, 'init')
    return;
end

updateall(hObject, eventdata, handles);



% --- Executes during object creation, after setting all properties.
function blankchannelmode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to blankchannelmode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in updatebutton.
function updatebutton_Callback(hObject, eventdata, handles)
% hObject    handle to updatebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateall(hObject, eventdata, handles);


% --- Executes on button press in loadbutton.
function loadbutton_Callback(hObject, eventdata, handles)
% hObject    handle to loadbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

aum = getappdata(handles.NexusChannels, 'AmpsUsedMode');

[selchanset, ChanNames, Montage] = nexus_prompt_select_montage ('nexuschannels');
if ~isempty(selchanset)
    if length(ChanNames) <= 34
        switch aum
            case 1
                NexusChanNames.nexus1 = ChanNames;
                NexusChanNames.nexus2 = Montage.numericnexus2;
            case 2
                NexusChanNames.nexus1 = ChanNames;
                NexusChanNames.nexus2 = Montage.numericnexus2;
            case 3
                NexusChanNames.nexus1 = Montage.numericnexus1;
                NexusChanNames.nexus2 = ChanNames;
        end
    else
        NexusChanNames.nexus1 = ChanNames(1:33);
        NexusChanNames.nexus2 = ChanNames(34:end);
    end
    loadallchans (NexusChanNames, handles);
    updateall(hObject, eventdata, handles);
end


% --- Executes when user attempts to close NexusChannels.
function NexusChannels_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to NexusChannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

updateall(hObject, eventdata, handles);
delete(hObject);


% --- Executes on button press in swapamps.
function swapamps_Callback(hObject, eventdata, handles)
% hObject    handle to swapamps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ch1gnd = get(handles.ch1gnd, 'String');
ch1chn = cell(1,32);
for chi = 1:32
    ch1chn{chi} = get(handles.(['ch1' num2str(chi, '%02i')]), 'String');
end

ch2gnd = get(handles.ch2gnd, 'String');
ch2chn = cell(1,32);
for chi = 1:32
    ch2chn{chi} = get(handles.(['ch2' num2str(chi, '%02i')]), 'String');
end

for chi = 1:32
    set(handles.(['ch1' num2str(chi, '%02i')]), 'String', ch2chn{chi});
end

for chi = 1:32
    set(handles.(['ch2' num2str(chi, '%02i')]), 'String', ch1chn{chi});
end
