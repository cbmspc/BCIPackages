% Classwise and Traditional (combined) Principal Component Analysis based
% on work done by Zoran Nenadic. Cut-off criteria: Half of mean of non-zero
% eigenvalues

function DRmatC = dataproc_func_ctpca(TrainData, TrainLabels, DRparm)

Nobs = length(TrainLabels);
if Nobs ~= size(TrainData,1)
    error('Number of observations from TrainData and TrainLabels disagree.');
end

classes = unique(TrainLabels);

Nclass = length(classes);

for c = 1:Nclass
    NtrialA(c) = length(find(TrainLabels == classes(c)));
end

for c = 1:Nclass
    ClassData = TrainData(find(TrainLabels==classes(c)),:);
    [coeff, latent] = dataproc_func_princomp(ClassData);
    coeffrC{c} = coeff(:,find( latent > mean(latent) ));
    sampmu{c} = mean(ClassData,1);
end

% Calculate traditional PCA
[coeff, latent] = dataproc_func_princomp(TrainData);
coeffrC{Nclass+1} = coeff(:,find( latent > mean(latent) ));

sampmuall = mean(TrainData,1);

% Calculate between-class covariance
for c = 1:Nclass
    Data_b(c,:) = sqrt(NtrialA(c) / Nobs) * (sampmu{c} - sampmuall);
end

W_b = dataproc_func_princomp(Data_b);

% Calculate principle subspace basis
for c = 1:length(coeffrC)
    DRmatC{c} = orth([coeffrC{c} W_b]);
end

