function [PPYoutput,f,TTYoutput] = fftplot(y, Fs, Frange, YScaling, Res)

BW = Fs/2;

if isempty(who('Frange')) || isempty(Frange)
    Frange = [0 BW];
end

% if isempty(who('styles')) || isempty(styles)
%     styles = '';
% end

if ~exist('Res','var') || isempty(Res)
    FFTFS = size(y,1);
    if mod(FFTFS,2)
        FFTFS = FFTFS+1;
    end
else
    FFTFS = 2^ceil(log2(BW)+1)*Res;
end

if ~exist('YScaling','var') || isempty(YScaling)
    YScaling = 'log';
end

NumFinalPoints = FFTFS/2+1;
Nc = size(y,2);
Ncolumn = min(2,Nc);
PlotsPerFig = 8;

% Consider possible user flipping
if size(y,1)/Fs < 0.3 || size(y,2) > 1024
    fprintf('It seems you have %i sample points (%.3g sec) and %i channels.\n', size(y,1), size(y,1)/Fs, size(y,2));
    fprintf('Did you flip the dimension? Signals should be in the first dimension.\n');
    fprintf('Ctrl+C now to stop operation. Otherwise, continuing without changing.\n');
    pause(10.0);
end

PPYoutput = zeros(NumFinalPoints,Nc);
TTYoutput = zeros(NumFinalPoints,Nc);

for c = 1:Nc
    Y = fft(y(:,c),FFTFS);
    PY = Y .* conj(Y) / FFTFS;
    TY = atan2(imag(Y),real(Y));

    f = Fs * (0:FFTFS/2) / FFTFS;

    PPY = PY(1:FFTFS/2+1);
    TTY = TY(1:FFTFS/2+1);
    PPYoutput(:,c) = PPY;
    TTYoutput(:,c) = TTY;

    if nargout == 0
        fid = intersect(find(f<=Frange(2)),find(f>=Frange(1)));
        fp = f(fid);
        PPYp = PPY(fid);
        TTYp = TTY(fid);
        fignum = floor((c-1)/PlotsPerFig)+1;
        figure(fignum);
        subplot(ceil(min(PlotsPerFig,Nc)/Ncolumn),Ncolumn,c-(fignum-1)*PlotsPerFig);
        switch YScaling
            case 'linear'
                plot(fp, PPYp);
            otherwise
                semilogy(fp,PPYp);
        end
        title(['PSD (' num2str(c) ')']);
        xlabel('Frequency (Hz)');
        ylabel('Power');
    end
end

if length(Frange) == 2
    fid = intersect(find(f<=Frange(2)),find(f>=Frange(1)));
    f = f(fid);
    PPYoutput = PPYoutput(fid,:);
    TTYoutput = TTYoutput(fid,:);
end
