% (traditional) Principal Component Analysis
% Cut-off criteria: Half of mean of non-zero eigenvalues or as specified in
% DRparm

function [DRmat, latent] = dataproc_func_pca(TrainData, TrainLabels, DRparm) %#ok<INUSL>
[coeff, latent] = dataproc_func_princomp(TrainData);

% if length(who('DRparm')) & length(DRparm) & isscalar(DRparm) & DRparm > 0
%     DRmat = coeff(:,1:DRparm);
% else
%     DRmat = coeff(:,find( latent > mean(latent) ));
% end
% 2011-11-30 more efficient coding
% if exist('DRparm','var') && ~isempty(DRparm) && isscalar(DRparm) && DRparm > 0
%     DRmat = coeff(:,1:DRparm);
% else
%     DRmat = coeff(:, latent > mean(latent) );
% end


% 2016-04-11 more options
if exist('DRparm','var') && ~isempty(DRparm) && isscalar(DRparm) && DRparm >= 0
    if DRparm == 0
        % Use the default: Keep above non-zero mean of eigenvalues
        DRmat = coeff(:, latent > mean(latent) );
    elseif DRparm < 1
        % Keep fraction of the sum of eigenvalues
        DRmat = coeff(:, 1:find(cumsum(latent)/sum(latent) >= DRparm,1) );
    else
        % Keep the number of eigenvectors
        DRmat = coeff(:, 1:min(length(latent),DRparm) );
    end
end