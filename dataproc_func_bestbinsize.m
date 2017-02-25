function Binsize = bestbinsize (Nsamp, Nchan, MinNtrial)

for b = 1:Nsamp+1
    Nbin = floor(Nsamp/b);
    if Nbin * Nchan < MinNtrial
        break;
    end
end

for c = b:Nsamp+1
    if floor(Nsamp/c) < Nbin
        c = c-1;
        break;
    end
end

Binsize = c;

if Binsize > Nsamp
    warning('Not enough number of trials to bin the data.');
end

% DataDim = Nsamp .* Nchan;
% 
% Nbin = floor(Nsamp ./ Binsize);
% Nlostsamp = Nsamp - Nbin .* Binsize;
% BinDim = Nbin .* Nchan;
% 
% disp(['    Nsamp = ' num2str(Nsamp)]);
% disp(['    Nchan = ' num2str(Nchan)]);
% disp(['MinNtrial = ' num2str(MinNtrial)]);
% disp(['  DataDim = ' num2str(DataDim)]);
% disp(['  Binsize = ' num2str(Binsize)]);
% disp(['     Nbin = ' num2str(Nbin)]);
% disp(['Nlostsamp = ' num2str(Nlostsamp)]);
% disp(['   BinDim = ' num2str(BinDim)]);
% 
