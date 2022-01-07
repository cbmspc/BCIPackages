function [eq, sub1equal, sub2equal] = subspace_equals(basis_set1, basis_set2, tol)
% basis_set1 is a set of tall vectors spanning a subspace
% basis_set2 is a set of tall vectors spanning a subspace
%
% If each vector inside basis_set1 can be written as a linear combination
% of basis_set2, and vice versa, the subspace spanned by basis_set1 is the
% same subspace spanned by basis_set2
%
% They don't need to be orthonormalized for this function to work.

if ~exist('tol','var') || isempty(tol)
    tol = eps * size(basis_set1,1);
end

sub1equal = false(1,size(basis_set1,2));
sub2equal = false(1,size(basis_set2,2));

for i = 1:size(basis_set1,2)
    % Can the ith vector in basis_set1 be written as v1 w1 + v2 w2 + ...
    % where v1 v2 ... are vectors of basis_set2 and w1 w2 ... are scalar
    % weights? The answer is yes if we can find w1 w2 ...
    % This is simply an Ax=b problem
    w = basis_set2 \ basis_set1(:,i);
    if norm(basis_set1(:,i) - basis_set2*w) < tol
        sub1equal(i) = true;
    end
end

for i = 1:size(basis_set2,2)
    % Can the ith vector in basis_set2 be written as v1 w1 + v2 w2 + ...
    % where v1 v2 ... are vectors of basis_set1 and w1 w2 ... are scalar
    % weights? The answer is yes if we can find w1 w2 ...
    % This is simply an Ax=b problem
    w = basis_set1 \ basis_set2(:,i);
    if norm(basis_set2(:,i) - basis_set1*w) < tol
        sub2equal(i) = true;
    end
end

eq = min([sub1equal sub2equal]) == 1;
    
