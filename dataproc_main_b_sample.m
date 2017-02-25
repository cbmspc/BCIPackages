%% Main program for data processing on sample (online data)
% function [Sdimdata_sample] = dataproc_main_training (Binsize, Smat,
% SourceData)
%
% Inputs:
% Binsize and Smat are obtained from calling dataproc_main_training on the
% training data. They must be the same for both training and online sample
% SourceData = Matrix of Windowed Raw Data (#chans, #samps) or 
% (#chans, #samps, #trials)
% Output: Sdimdata_sample is ready to be transformed by Fmat.
%

function [Sdimdata_sample] = dataproc_main_testing (Binsize, Smat, SourceData)

%% Gather the dimensions
if iscell(SourceData) || length(size(SourceData)) < 2
    error('SourceData must be a 2D or 3D matrix');
end

%% Binning and reshaping
DimData = dataproc_func_binning(SourceData,Binsize);

%% Sphereing transformation
Sdimdata_sample = Smat * DimData;
clear('DimData');

%% Example code
% Fs = 200;
% M = 2; % Feature space dimension
% [Fmat, Flabel, Smat, Sinvmat, Binsize, Sdimdata_training] = dataproc_main_training(Fs, M, {TrainingEEG.moveleft, TrainingEEG.moveright});
% Ftrain = Fmat * Sdimdata_training;
% Sdimdata_online = dataproc_main_sample(Binsize, Smat, OnlineEEG);
% Fonline = Fmat * Sdimdata_online;
% CQ = classify(Fonline', Ftrain', Flabel, 'quadratic', 'empirical');
% switch CQ
%     case 0
%         disp('decoded to be LEFT');
%     case 1
%         disp('decoded to be RIGHT');
% end
