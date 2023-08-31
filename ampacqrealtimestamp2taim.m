% from AmpAcqRealTimestamp to datenum
% TAIm: Approximate TAIM of the acquisition. Only need the date component. Can be off by up to a week.
% Example
%{
  tmin = 34490749
  tms = 35517
  aarts = bitand(tms,65535) + bitshift(bitand(floor(mod(tmin,1440)),2^12-1),16) + bitshift(bitand(mod(floor(tmin/1440),32),2^6-1),27)
  taim_expmdate = datenum2taim(datenum('2023-07-30')); % Enter the experiment date
  [tmin_expm, tms_expm, tstr_expm, datenum_expm, datestr_expm] = ampacqrealtimestamp2taim(aarts, taim_expmdate)
%}
% Reminder that TAI is in zulu time
function [TAIminute, TAImillisec, TAImmString, DateNum, ISODateStr] = ampacqrealtimestamp2taim(AmpAcqRealTimestamp, TAIm)
% LSB 16 bits: milliseconds
% Next 11 bits: minutes in the day
% Last 5 bits: days % 32

TAImillisec = bitand(AmpAcqRealTimestamp, 65535);
minutes_in_day = bitshift(bitand(AmpAcqRealTimestamp, 134152192), -16);
days_mod_32 = bitshift(bitand(AmpAcqRealTimestamp, 4160749568), -27);
TAIm = floor(TAIm);
TAId = floor(TAIm / 1440);
m = mod(TAId, 32);

if m == days_mod_32
    TAIminute = TAId * 1440 + minutes_in_day;
elseif mod(m + 1,32) == days_mod_32
    TAIminute = TAId * 1440 + minutes_in_day + 1440;
elseif mod(m - 1,32) == days_mod_32
    TAIminute = TAId * 1440 + minutes_in_day - 1440;
elseif mod(m + 2,32) == days_mod_32
    TAIminute = TAId * 1440 + minutes_in_day + 1440*2;
elseif mod(m - 2,32) == days_mod_32
    TAIminute = TAId * 1440 + minutes_in_day - 1440*2;
elseif mod(m + 3,32) == days_mod_32
    TAIminute = TAId * 1440 + minutes_in_day + 1440*3;
elseif mod(m - 3,32) == days_mod_32
    TAIminute = TAId * 1440 + minutes_in_day - 1440*3;
elseif mod(m + 4,32) == days_mod_32
    TAIminute = TAId * 1440 + minutes_in_day + 1440*4;
elseif mod(m - 4,32) == days_mod_32
    TAIminute = TAId * 1440 + minutes_in_day - 1440*4;
elseif mod(m + 5,32) == days_mod_32
    TAIminute = TAId * 1440 + minutes_in_day + 1440*5;
elseif mod(m - 5,32) == days_mod_32
    TAIminute = TAId * 1440 + minutes_in_day - 1440*5;
elseif mod(m + 6,32) == days_mod_32
    TAIminute = TAId * 1440 + minutes_in_day + 1440*6;
elseif mod(m - 6,32) == days_mod_32
    TAIminute = TAId * 1440 + minutes_in_day - 1440*6;
elseif mod(m + 7,32) == days_mod_32
    TAIminute = TAId * 1440 + minutes_in_day + 1440*7;
elseif mod(m - 7,32) == days_mod_32
    TAIminute = TAId * 1440 + minutes_in_day - 1440*7;
else
    error('Provided TAIM is too far off from the acquisition date');
end

TAImmString = sprintf('%d:%05d', TAIminute, TAImillisec);

DateNum = taim2datenum(TAImmString);
ISODateStr = [datestr(DateNum, 'yyyy-mm-dd HH:MM:SS.fff') ' UTC' num2str(tzoffsethour())] ;

