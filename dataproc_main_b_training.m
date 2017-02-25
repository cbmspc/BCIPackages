%% Main program for data processing on training
% function [Fmat, Flabel, Smat, Sinvmat, Binsize, Sdimdata] =
% dataproc_main_training (SamplingFrequency, Fdim, CellSourceData, OverrideBinsize)
%
% Input: 1D Cells of Windowed Raw Data {(#chans, #samps, #trials), (#chans,
% #samps, #trials), ...} each element of cell represent one class
%
% Fdim = Featurespace dimension
%
% OverrideBinsize (OPTIONAL): If specified, will override the binsize
% instead of using automatic determination.
%
% Output: Feature extraction matrix, labels (labels are integers starting
% from 0), sphering matrix, desphering matrix, bin size used, and the sphered
% reshaped binned data
%
% The end of this file contains instructions to process online data



function [Fmat, Flabel, Smat, Sinvmat, Binsize, Sdimdata] = dataproc_main_training (CellSourceData, Fdim, OverrideBinsize)

%% Pre-populate output variables
Fmat = [];
Ftrain = [];
Flabel = [];
Smat = [];
Sinvmat = [];
Binsize = [];

%% Gather the dimensions
if ~iscell(CellSourceData)
    error('CellSourceData must be a cell');
end
Nclass = length(CellSourceData);
if Nclass == 0
    error('Size of CellSourceData is zero. Nothing to process!');
end
Nchan = size(CellSourceData{1},1);
Nsamp = size(CellSourceData{1},2);
for c = 1:Nclass
    VNchan(c) = size(CellSourceData{c},1);
    VNsamp(c) = size(CellSourceData{c},2);
    VNtrial(c) = size(CellSourceData{c},3);
    % Create labels
    Flabel = cat(2,Flabel,(c-1)*ones(1,VNtrial(c)));
end
% Make sure #chans and #samps are the same
if max(abs(diff(VNchan)))
    error('Number of channels is different among classes');
end
if max(abs(diff(VNsamp)))
    error('Number of samples is different among classes');
end
clear('VNchan','VNsamp');

%% Determine the best bin size. 
% disp('Beginning binsize determination..');
% Binned data dimension must be less than the lowest number of trials among
% all classes so that the covariance matrix stays full rank.

% Determine minimum number of trials among classes. This will be the
% limiting factor for bin size (less trials = bigger bin size)
MinNtrial = min(VNtrial);
% DataDim = Nchan * Nsamp;
% for Binsize = 1:Nsamp
%     % If binsize does not divide WindowDim, skip
%     % If WindowDim is prime number, don't bother
%     if ~isprime(Nsamp) & mod(Nsamp, Binsize)
%         continue;
%     end
%     % For covariance matrix to be full rank, the Binned data dimension must
%     % be less than the lowest number of trials among all classes
%     if DataDim / Binsize < MinNtrial-1
%         break;
%     end
% end
% 
% if Binsize == Nsamp
%     warning('Problem determining the best Binsize');
% end

if nargin >= 3
    % Override Binsize with this value
    Binsize = OverrideBinsize;
else
    % Use the smallest possible Binsize while minimizing loss of time
    % points at the end of window.
    Binsize = dataproc_func_bestbinsize(Nsamp, Nchan, MinNtrial);
end

if Binsize > Nsamp || Binsize <= 0
    error('Not enough number of trials to create a full rank covariance matrix: Reduce number of channels or size of time window.');
end

%% Binning, reshaping, and concatenate all classes together along trials
% disp('Beginning data binning, reshaping, and concatenating..');
DimData = [];
for c=1:Nclass
    DimData = cat(2,DimData,dataproc_func_binning(CellSourceData{c},Binsize));
end

%% Sphereing transformation
% disp('Beginning Sphering transformation..');
[Sdimdata, Smat, Sinvmat] = dataproc_func_sphereing(DimData);
clear('DimData');

%% Feature extraction
% disp('Beginning Feature extraction..');
Fmat = dataproc_func_featureextraction(Fdim, Sdimdata, Flabel);

%% Instructions to process online data
% 1. SourceData = Window and preprocess data as necessary
% 2. DimData = bin, reshape, and concatenate
% 3. Sdimdata_online = sphereing transformation
% 4. Fonline = Fmat * Sdimdata_online
% 5. Ftrain = Fmat * Sdimdata
% 6. Use your favorite classifier using Fonline and Ftrain

