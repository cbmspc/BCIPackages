function [newvectors, sphereingMatrix, desphereingMatrix] = sphereing(vectors)
% Performs sphereing transformation and returns transormed data vectors
% along with transformation matrices.
%
% vectors - data in row vectors {m x n}
% sphereingMatrix - sphereing transformation matrix {mxm}
% desphereingMatrix - reverse transformation matrix {mxm}
%
covarianceMatrix = cov(vectors', 1);

% Calculate the eigenvalues and eigenvectors of covariance
% matrix.
[E, D] = eig(covarianceMatrix);

% ========================================================
% Calculate the whitening and dewhitening matrices (these handle
% dimensionality simultaneously).
sphereingMatrix = inv(sqrt(D)) * E';
desphereingMatrix = E * sqrt(D);

newvectors = sphereingMatrix * vectors;

%cov(newvectors',1)   % $$$ test - must be diagonal with ones

return;