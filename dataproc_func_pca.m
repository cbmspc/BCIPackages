% Dimension Reduction using traditional Principal Component Analysis
%
% DRparm = 0 : (Default) Half of mean of non-zero eigenvalues
%          Negative fraction less than 1 : Retain eigenvalues with
%              magnitude at least this fraction of largest eigenvalue
%          Positive fraction less than 1 : Fraction of energy retained
%          Positive integer >= 1 : Number of retained eigenvectors
%
% Output:
%   DRmat = Dimension reduction transformation matrix, i.e. the matrix
%           formed by the retained Principal Axes
%  latent = Eigenvalues of the corrsponding principal axes


function [DRmat, latent] = dataproc_func_pca(TrainData, TrainLabels, DRparm) %#ok<INUSL>
if any(isnan(TrainData(:)))
    error('TrainData contains NaN');
end

[coeff, latent] = dataproc_func_princomp(TrainData);


if ~exist('DRparm','var') || isempty(DRparm) || ~isscalar(DRparm)
    DRparm = 0;
end

if DRparm == 0
    % Use the default: Keep above non-zero mean of eigenvalues
    DRmat = coeff(:, latent > mean(latent(latent>0)) );
elseif DRparm < 0 && DRparm > -1
    DRmat = coeff(:, (latent/max(latent) > -DRparm) );
elseif DRparm <= -1
    % Undefined behavior. Use default
    DRmat = coeff(:, latent > mean(latent(latent>0)) );
elseif DRparm < 1
    % Keep fraction of the sum of eigenvalues
    DRmat = coeff(:, 1:find(cumsum(latent)/sum(latent) >= DRparm,1) );
else
    % Keep the number of eigenvectors
    % This number cannot exceed the number of non-zero eigenvalues
    DRmat = coeff(:, 1:min(nnz(latent>0),DRparm) );
end
