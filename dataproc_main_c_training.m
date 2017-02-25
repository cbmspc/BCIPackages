% function [FPmatC, DimDataC, PCAbasis, FmatC] = dataproc_main_c_training
% (CellSourceData, Fdim)
%
% Training program based on classwise principle component analysis (CPCA)
% and approximated information discriminant analysis (AIDA).
%
% FPmatC is the combination of PCA basis matrix and IDA feature matrix.
%
% The counterpart of this program is dataproc_main_c_classify for
% classifying test or unknown data

function [FPmatC, DimDataC, PCAbasis, FmatC] = dataproc_main_c_training (CellSourceData, Fdim)

Flabel = [];

%% Gather the dimensions
if ~iscell(CellSourceData)
    error('CellSourceData must be a cell');
end
Nclass = length(CellSourceData);
if Nclass == 0
    error('Size of CellSourceData is zero. Nothing to process!');
end
% Nchan = size(CellSourceData{1},1);
% Nsamp = size(CellSourceData{1},2);
Vreshapedone = zeros(1,Nclass);
for c = 1:Nclass
    if length(size(CellSourceData{c})) == 3
        VNchan(c) = size(CellSourceData{c},1);
        VNsamp(c) = size(CellSourceData{c},2);
        VNtrial(c) = size(CellSourceData{c},3);
    else  % Data is already reshaped
        Vreshapedone(c) = 1;
        VNchan(c) = 1;
        VNsamp(c) = size(CellSourceData{c},1);
        VNtrial(c) = size(CellSourceData{c},2);
    end
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

%% Reshaping, Principal Component Analysis, and PCA reduction
DimDataC = [];
DimDataA = [];
for c = 1:Nclass
    if Vreshapedone(c)    % Already reshaped
        DimDataC{c} = CellSourceData{c};
    else                  % Reshape
        for t = 1:size(CellSourceData{c},3)
            DimDataC{c}(:,t) = reshape(CellSourceData{c}(:,:,t), 1, []);
        end
    end
    DimDataA = cat(2,DimDataA,DimDataC{c});

    % PCA on each class
    [coeff, latent] = dataproc_func_princomp(DimDataC{c}.');
    coeffrC{c} = coeff(:,find( latent > mean(latent) ));

    sampmu{c} = mean(DimDataC{c},2);
    %size(coeffrC{c})
end

% % PCA on all classes
% [coeff, latent] = dataproc_func_princomp(DimDataA');
% coeffrC{Nclass+1} = coeff(:,find( latent > mean(latent) ));


clear('score','latent','coeff','Vreshapedone');

sampmuall = mean(DimDataA,2);
% Sigma_b = (sampmu{1} - sampmuall) * (sampmu{1} - sampmuall)' * VNtrial(1) / sum(VNtrial);
% for c = 2:Nclass
%     Sigma_b = Sigma_b + (sampmu{c} - sampmuall) * (sampmu{c} - sampmuall)' * VNtrial(c) / sum(VNtrial);
% end
% W_b = orth(Sigma_b);

% Use intelligent PCA to calculate eigenvectors of between-class covariance
NtrialTotal = sum(VNtrial);
for i = 1:Nclass
    Data_b(i,:) = sqrt(VNtrial(i) / NtrialTotal) * (sampmu{i} - sampmuall);
end

W_b = dataproc_func_princomp(Data_b);


for c = 1:Nclass
    
    % Calculate principle subspace basis
    PCAbasis{c} = orth([coeffrC{c} W_b]);
    
    % Principle subspace projection of all observations
    % There are Nclass subspaces that the same data can project onto
    DimData_PA_Proj{c} = ( DimDataA.' * PCAbasis{c} ).';
    
    % Feature extraction on each classwise subspace
    FmatC{c} = dataproc_func_featureextraction(Fdim, DimData_PA_Proj{c}, Flabel, 0);
    FPmatC{c} = FmatC{c} * PCAbasis{c}.';
end

% % PCA on all classes
% PCAbasis{Nclass+1} = orth([coeffrAll W_b]);
% DimData_PA_Proj{Nclass+1} = ( DimDataA.' * PCAbasisAll ) .';
% FmatC{Nclass+1} = dataproc_func_featureextraction(Fdim, DimData_PA_ProjAll, Flabel);
% FPmatC{Nclass+1} = FmatC{Nclass+1} * PCAbasis{Nclass+1}.';


clear('Sigma_b','W_b','coeffrC');


