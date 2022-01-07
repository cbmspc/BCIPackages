% Create Biquadratic coefficients

function Hd = bqc_create(Fs, Fc1, Fc2, N, SWsilent)

if ~exist('SWsilent','var') || isempty(SWsilent)
    SWsilent = 0;
end

%Fs = 500;  % Sampling Frequency

%N   = 4;   % Order
%Fc1 = 8;   % First Cutoff Frequency
%Fc2 = 12;  % Second Cutoff Frequency

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.bandpass('N,F3dB1,F3dB2', N, Fc1, Fc2, Fs);
Hd = design(h, 'butter');

% Modified from https://tomlankhorst.nl/filter-controller-cpp-implementation/
sos = Hd.sosMatrix;
scales = Hd.scaleValues;
overall_scale = prod(Hd.scaleValues);

% for i = 1:size(sos,1)
%     sos(i,[1:3]) = sos(i,[1:3]) * scale(i);
% end


% i = 0;
% for s = sos.'
%     i = i + 1;
%     fprintf('BiQuad bq%d( %.16g, %.16g, %.16g, %.16g, %.16g );\n', i, s(1), s(2), s(3), s(5), s(6));
% end

if ~SWsilent
    i = 0;
    in = size(sos,1);
    for s = sos.'
        i = i + 1;
        fprintf('{ %.16g, %.16g, %.16g, %i, %i, %.16g, %.16g, %.16g, %.16g, %.16g, %.16g };\n', Fs, Fc1, Fc2, in, i-1, s(5), s(6), scales(i), s(1), s(2), s(3));
    end

    fprintf('#define FILTER_SOS_BP_NUMER { %.16g, %.16g, %.16g }\n', sos(1,1), sos(1,2), sos(1,3));
    fprintf('#define FILTER_SOS_%iHZ_%i_%i_BP_ORD%i_DENOM { ', Fs, Fc1, Fc2, N);
    for i = 1:size(sos,1)
        fprintf('%.16g, %.16g, ', sos(i,5), sos(i,6));
    end
    fprintf('\b\b }\n');
    fprintf('#define FILTER_SOS_%iHZ_%i_%i_BP_ORD%i_SCALE { ', Fs, Fc1, Fc2, N);
    for i = 1:size(sos,1)
        fprintf('%.16g, %.16g, ', scales(i));
    end
    fprintf('\b\b }\n');
    %fprintf('#define FILTER_SCALE_%iHZ_%i_%i_BP_ORD%i %.16g\n', Fs, Fc1, Fc2, N, overall_scale);
    fprintf('\n\n');
end



% fprintf('\nbqc');
% i = 0;
% for s = sos.'
%     i = i + 1;
%     fprintf('.add( &bq%d )', i);
% end
% fprintf(';\n');