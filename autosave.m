% Sourced from c:\bin\matlab_packages\autosave.m
function tout = autosave(varargin)
% AUTOSAVE automatically saves data in the base workspace
% AUTOSAVE saves all of the variables in the base workspace every 10
% minutes, to the matlab.mat file.
% AUTOSAVE(PERIOD) will save data every PERIOD minutes to the matlab.mat file.
% AUTOSAVE(PERIOD,SAVEFILE) will save variables to the file whose name is
%  specified by SAVEFILE.
% AUTOSAVE(PERIOD,SAVEFILE,SAVEARGS) will pass the options in the cell
%  array as inputs to the SAVE function.
% AUTOSAVE STOP will disable autosaving.
% AUTOSAVE START will restart autosaving.
% AUTOSAVE DELETE will disable autosaving, delete the autosave file, and
%  remove the timer used to perform the saving.
%
% Examples:
%
%    % autosave with defaults
%
%    autosave
% 
%
%    % autosave every 3 minutes
%
%    autosave(3)
%
%
%    % autosave every 7 minutes, to data.mat
%
%    autosave(7,'data.mat')
%
%
%    % save variables whose names start with 'x' into 'test.mat'
%    % every 5 minutes.
%
%    autosave(5,'test.mat','x*')
% 
%{
Copyright (c) 2016, The MathWorks, Inc.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the distribution.
    * In all cases, the software is, and all modifications and derivatives
      of the software shall be, licensed to you solely for use in conjunction
      with MathWorks products and service offerings.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
%}



% defaults
period = 30;
filename = [feature('logdir') filesep 'matlab.mat'];
saveargs = {'-v7.3', '-nocompression'};
minimum_freespace_fraction_required = 0.25;

% check inputs
error(nargchk(0,3,nargin));

% handle STOP/START/DELETE
if nargin == 1 && ischar(varargin{1})
    switch  lower(varargin{1})
        case 'start'
            t = getappdata(0,'AutoSaveTimerHandle');
            if isempty(t)
                t = timerfind('Tag','AutoSaveTimer');
            end
            if isempty(t) || ~isvalid(t)
                t = autosave;
            end
            if strcmp(get(t,'Running'),'off')
                start(t);
            end
        case 'stop'
            t = getappdata(0,'AutoSaveTimerHandle');
            if isempty(t)
                t = timerfind('Tag','AutoSaveTimer');
            end
            if ~isempty(t) && isvalid(t)
                stop(t);
            end            
        case 'delete'
            t = autosave('stop');
            if isvalid(t)
                delete(t);
            end
            filename = getappdata(0,'AutoSaveFile');
            if ~isempty(filename) && ischar(filename) && exist(filename, 'file')
                delete(filename);
            end            
    end
    
    % return early
    if nargout
        tout = t;
    end
    return;
end


% assign inputs
switch nargin
    case 1
        period = varargin{1};
    case 2
        period = varargin{1};
        filename = varargin{2};
    case 3 
        period = varargin{1};
        filename = varargin{2};
        saveargs = varargin{3};
end

% check inputs. 
msg = checkPeriod(period); error(msg);
msg = checkFilename(filename); error(msg);
msg = checkSaveArgs(saveargs); error(msg);

% convert to seconds
periodSeconds = period * 60;

t = getappdata(0,'AutoSaveTimerHandle');

if isempty(t) || ~isvalid(t)
    % create timer object
    t = timer('TimerFcn',{@savedata,filename,minimum_freespace_fraction_required,saveargs},...
        'Period',periodSeconds,...
        'ExecutionMode','fixedSpacing',...
        'Tag','AutoSaveTimer');
    % save timer handle
    setappdata(0,'AutoSaveTimerHandle',t);
    setappdata(0,'AutoSaveFile',filename);
end

      
% start the timer
if strcmp(get(t,'Running'),'off')
    start(t);
end

if nargout
    tout = t;
end

%%%%%%%%%%%  save data  %%%%%%%%%%%

function savedata(h,ev,filename,minimum_freespace_fraction_required,saveargs)
%SAVEDATA saves data from the base workspace

% convert cell array to a single string
if iscellstr(saveargs)
    tmpsaveargs = '';
    for n = 1:length(saveargs)
        tmpsaveargs = [tmpsaveargs ' ' saveargs{n}];
    end
    saveargs = tmpsaveargs;
end

% For matlabrc autosave, automatically delete >1.5-day old autosaves after MATLAB uptime > 2.5 hours
[Folder, Filenamepart] = fileparts(filename);
if isequal(Folder,feature('logdir')) && ~isempty(regexp(Filenamepart, '^matlabbaseautosave\d+', 'match', 'once'))
    if cputime > 9000
        tmp_list = dir([Folder filesep regexprep(Filenamepart, '\d+', '*')]);
        for i = 1:length(tmp_list)
            if now - tmp_list(i).datenum > 1.5 && ~isempty(regexp(tmp_list(i).name, '^matlabbaseautosave\d+\.mat$', 'match', 'once'))
                delete([Folder filesep tmp_list(i).name]);
            end
        end
    end

    % while handling matlabrc autosave, calculate space requirement
    tmp_whos = evalin('base', 'whos');
    tmp_bytes = sum([tmp_whos.bytes]);
    tmp_req = tmp_bytes / minimum_freespace_fraction_required;
    FileObj = java.io.File(Folder);
    usable_bytes = FileObj.getUsableSpace;

    if tmp_req > usable_bytes
        % Try to more aggressively delete old autosaves
        tmp_list = dir([Folder filesep regexprep(Filenamepart, '\d+', '*')]);
        for i = 1:length(tmp_list)
            if now - tmp_list(i).datenum > 2/24 && ~isempty(regexp(tmp_list(i).name, '^matlabbaseautosave\d+\.mat$', 'match', 'once'))
                delete([Folder filesep tmp_list(i).name]);
                FileObj = java.io.File(Folder);
                usable_bytes = FileObj.getUsableSpace;
                if tmp_req <= usable_bytes
                    break;
                end
            end
        end
    end
    if tmp_bytes >= usable_bytes
        % Cannot autosave after all.
        return
    end


end

% create string to execute that will save the data
savestr = ['try save ' filename ' ' saveargs ';end'];
evalin('base',savestr);




%%%%%%%%%%  check inputs  %%%%%%%%%%

function msg = checkPeriod(period)
msg = '';
if ~isnumeric(period) || numel(period) ~= 1
    msg = 'First input must be a numeric scalar.';
end

function msg = checkFilename(filename)
msg = '';
if ~ischar(filename) || size(filename,1) ~= 1
    msg = 'Second input must be a string array.';
end

function msg = checkSaveArgs(saveargs)
msg = '';
if isempty(saveargs)
    return
end
if ~ischar(saveargs) && ~iscellstr(saveargs)
    msg = 'Third input must be a string array or cell array of strings.';
end