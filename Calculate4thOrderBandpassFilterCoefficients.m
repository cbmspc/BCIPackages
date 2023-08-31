function [a1, a2, scale] = Calculate4thOrderBandpassFilterCoefficients(Fs, fL, fH)
% [~, a1, a2, scale] = bqc_create(Fs, Fc1, Fc2, 4, 1); % 4th order (meant for 2 iterations)

%Derived by Dr. Zoran Nenadic on 2023-08-25
%Verified by Dr. Po T Wang on 2023-08-28. Results identical to Matlab

T = 1/Fs;

%pre-warp the corner frequencies
wL = 2*Fs*tan(fL*(2*pi/Fs)/2);
wH = 2*Fs*tan(fH*(2*pi/Fs)/2);

%introduce the (geometric) mean frequency and bandwidth 
wm = sqrt(wL*wH);
W = wH-wL;

%define the gain factor for the digital filter
K = 2/T;

% %Discretized 2nd filter (bilinear transform)
% num = [W*K/(K^2+W*K+wm^2), 0, -W*K/(K^2+W*K+wm^2)];
% den = [1, (2*wm^2-2*K^2)/(K^2+W*K+wm^2), (K^2-W*K+wm^2)/(K^2+W*K+wm^2)];

%test the 4th order Butterworth digital filter (Matlab) against the one 
% derived by the bilinear transform 

%alpha (alp), beta (bet), gamma (gam) and delta (del) necessary for the
%factorization of the 4th BP filter into a product of two biquads.
A = W/wm;
D = sqrt(A^4+16)+A^2;   %from Wolfram alpha
alp = 4*sqrt(2)*W/(D+A*sqrt(2*D)+4);
bet = 4*wm^2/(D+A*sqrt(2*D));
gam = sqrt(2)*W*(D+A*sqrt(2*D))/(D+A*sqrt(2*D)+4);
del = (D+A*sqrt(2*D))*wm^2/4;

%Discretized filters (bilinear transform) biquad #1 and #2

scale = [W*K/(K^2+alp*K+bet); 
    W*K/(K^2+gam*K+del)];

a1 = [(2*bet-2*K^2)/(K^2+alp*K+bet);
    (2*del-2*K^2)/(K^2+gam*K+del)];

a2 = [(K^2-alp*K+bet)/(K^2+alp*K+bet);
    (K^2-gam*K+del)/(K^2+gam*K+del)];

