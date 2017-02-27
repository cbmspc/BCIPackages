function XYoffsets = eeg_electrode_to_gridcoord (ChanNames)
% Convert from electrode name(s) to a list square grid coordinates
if isstr(ChanNames)
    ChanNames = string_to_cell(ChanNames,' ,-_');
end
Nchan = length(ChanNames);

XYoffsets = zeros(Nchan,2);

for i = 1:Nchan
    [XYoffsets(i,2),XYoffsets(i,1)] = e_convert(ChanNames{i});
end



function [Voffset, Hoffset] = e_convert (EName)
EName = upper(EName);

Hoffset = 0;  % + is right, - is left
Voffset = 0;  % + is anterior, - is posterior

% Interprete the horizontal offset
if strcmp(EName,'A1')
    Hoffset = -6;
elseif strcmp(EName,'M1')
    Hoffset = -6;
elseif strcmp(EName,'A2')
    Hoffset = 6;
elseif strcmp(EName,'M2')
    Hoffset = 6;
else
    switch EName(end)
        case 'Z'
            Hoffset = 0;
        case '1'
            Hoffset = -1;
        case '2'
            Hoffset = 1;
        case '3'
            Hoffset = -2;
        case '4'
            Hoffset = 2;
        case '5'
            Hoffset = -3;
        case '6'
            Hoffset = 3;
        case '7'
            Hoffset = -4;
        case '8'
            Hoffset = 4;
        case '9'
            Hoffset = -5;
        case '0'
            Hoffset = 5;
    end
end

% Interprete vertical offset
switch EName(1)
    case 'N'
        Voffset = 5;
    case 'F'
        switch EName(2)
            case 'P'
                Voffset = 4;
            case 'C'
                Voffset = 1;
            case 'T'
                Voffset = 1;
            otherwise
                Voffset = 2;
        end
    case 'A'
        switch EName(2)
            case 'F'
                Voffset = 3;
            otherwise
                Voffset = 0;
        end
    case 'C'
        switch EName(2)
            case 'P'
                Voffset = -1;
            otherwise
                Voffset = 0;
        end
    case 'T'
        switch EName(2)
            case 'P'
                Voffset = -1;
            otherwise
                Voffset = 0;
        end
    case 'P'
        switch EName(2)
            case 'O'
                Voffset = -3;
            otherwise
                Voffset = -2;
        end
    case 'O'
        Voffset = -4;
    case 'I'
        Voffset = -5;
    case 'M'
        Voffset = -1;
end


