% Dimension Reduction using Fukunaga-Koontz Transform
% This can only work on data with two classes
%
% Note: Only the first two classes (sorted in ascending order) are used
%
% DRparm is an array: [K1, K2, PCADISABLE]
%   K1 and K2 are positive integer or fraction. 
%     If fraction less than 1:
%        Keep EigVal greater than 1-K1 and less than K2
%     If integer greater or equal to 1:
%        Keep the first K1 EigVec and last K2 EigVec
%   PCADISABLE = 1 to disable using PCA to resolve rank deficient problem
%   (if any). Note: PCA is only used to condition data, not to discard any.
%   


function [DRmat, latent] = dataproc_func_fktdr(TrainData, TrainLabels, DRparm)

if (~exist('DRparm', 'var') || isempty(DRparm)), DRparm = [0 0 0]; end
DRparm = DRparm(:);
if (length(DRparm) == 1), DRparm(2) = DRparm(1); end
if (length(DRparm) == 2), DRparm(3) = 0; end
if (DRparm(1) < 0), DRparm(1) = 0; end
if (DRparm(2) < 0), DRparm(2) = 0; end
if (DRparm(3)), PCADISABLE = 1; else PCADISABLE = 0; end


classes = unique(TrainLabels);
N1 = nnz(TrainLabels == classes(1));
N2 = nnz(TrainLabels == classes(2));
Ndim = size(TrainData,2);
X1 = TrainData(TrainLabels == classes(1),:);
X2 = TrainData(TrainLabels == classes(2),:);


% FKT does not work when matrix is rank-deficient (e.g. channels with zero
% data). In this case, PCA is automatically performed on the data if not
% disallowed.
MATRANK = [ rank(X1,max(svd(X1))*1e-9), rank(X2,max(svd(X2))*1e-9) ];
MAXDIM = min(MATRANK);
if MAXDIM <= min([Ndim, N1, N2]) && ~PCADISABLE
    DRMPCA = dataproc_func_pca([X1;X2], [], MAXDIM);
else
    DRMPCA = eye(size([X1;X2],2));
end

X1 = X1*DRMPCA;
X2 = X2*DRMPCA;

S1=(X1'*X1)/N1;
S2=(X2'*X2)/N2;
[Phi,D1]=eig(S1+S2);
%W=Phi*sqrt(inv(D1));

% Better way to whiten that works for rank-deficient data:
tol = 1e-9;
W = Phi(:,diag(D1)>tol) * sqrt(inv(D1(diag(D1)>tol,diag(D1)>tol))) * Phi(:,diag(D1)>tol)';

% For X1
[Psi,D]=eig(W'*S1*W);
T=W*Psi;
latent=diag(D);
[latent, ix] = sort(latent, 'descend');
T = T(:,ix);


% Ideally, latent = 1-latent2.

keep = false(size(latent));
if DRparm(1) == 0 && DRparm(2) == 0
    % Use the default: Keep eigenvalues greater than 0.9 and less than 0.1
    %DRmat = T(:, [1:ceil(MATRANK(1)*0.10), end-ceil(MATRANK(2)*0.10)+1:end] );
    keep = keep | (latent>0.9 | latent<0.1);
else
    % Keep eigenvalues greater than 1-DRparm(1) and less than DRparm(2)
    % The smaller the threshold, the less dimensions are retained
    if DRparm(1) < 1
        keep = keep | (latent>(1-DRparm(1)));
    end
    if DRparm(2) < 1
        keep = keep | (latent<(DRparm(2)));
    end
    
    % Keep the number of eigenvectors from each end
    if DRparm(1) >= 1
        keep(1:DRparm(1)) = 1;
    end
    if DRparm(2) >= 1
        keep(end-DRparm(2)+1:end) = 1;
    end
end

DRmat = T(:,keep);
DRmat = DRMPCA*DRmat;


