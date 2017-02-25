% Classwise Principal Component Analysis based on work done by Zoran
% Nenadic. Cut-off criteria: Half of mean of non-zero eigenvalues or up to
% DRparm eigenvectors (if specified). Note: Final dimension size is usually
% DRparm plus up to the number of classes due to inclusion of between class
% covariances

function DRmatC = dataproc_func_cpca(TrainData, TrainLabels, DRparm)

Nobs = length(TrainLabels);
if Nobs ~= size(TrainData,1)
    error('Number of observations from TrainData and TrainLabels disagree.');
end

classes = unique(TrainLabels);

Nclass = length(classes);
Ndim = size(TrainData,2);

NtrialA = zeros(1,Nclass);
for c = 1:Nclass
    NtrialA(c) = length(find(TrainLabels == classes(c)));
end

if exist('DRparm','var') && ~isempty(DRparm) && isscalar(DRparm) && isnumeric(DRparm) && DRparm >= 0
else
    DRparm = 0;
end

coeffrC = cell(1,Nclass);
sampmu = cell(1,Nclass);
for c = 1:Nclass
    idc = find(TrainLabels==classes(c));
    [coeff, latent] = dataproc_func_princomp(TrainData(idc,:));
    if DRparm == 0
        % Use the default: Keep above non-zero mean of eigenvalues
        coeffrC{c} = coeff(:, latent > mean(latent) );
    elseif DRparm < 1
        % Keep fraction of the sum of eigenvalues
        coeffrC{c} = coeff(:, 1:find(cumsum(latent)/sum(latent) >= DRparm,1) );
    else
        % Keep the number of eigenvectors
        coeffrC{c} = coeff(:, 1:min(length(latent),DRparm) );
    end
    
    sampmu{c} = mean(TrainData(idc,:),1);
end

sampmuall = mean(TrainData,1);

% Calculate between-class covariance
Data_b = zeros(Nclass,Ndim);
for c = 1:Nclass
    Data_b(c,:) = sqrt(NtrialA(c) / Nobs) * (sampmu{c} - sampmuall);
end

W_b = dataproc_func_princomp(Data_b);

% Calculate principle subspace basis

DRmatC = cell(1,Nclass);
for c = 1:Nclass
    try
        DRmatC{c} = orth([coeffrC{c} W_b]);
    catch %#ok<CTCH>
        DRmatC{c} = zeros(size(coeffrC{c},1),0);
    end
end
