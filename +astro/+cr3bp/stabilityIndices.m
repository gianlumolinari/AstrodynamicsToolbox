function stab = stabilityIndices(M)
%STABILITYINDICES Compute eigenvalue-based stability diagnostics from monodromy matrix.
%
% INPUT
%   M : monodromy matrix [6x6]
%
% OUTPUT
%   stab : struct with fields
%       .eigvals          : eigenvalues of M
%       .eigvecs          : eigenvectors of M
%       .eigAbsSorted     : abs(eigenvalues), descending
%       .lambdaDominant   : dominant eigenvalue magnitude
%       .lambdaMin        : smallest eigenvalue magnitude
%       .nu               : simple planar-style stability index
%       .spectralRadius   : max(abs(eigvals))
%       .isApproximatelyStable : true if spectral radius ~ 1
%
% NOTES
%   For CR3BP periodic orbits, the full interpretation depends on family and
%   dimension. The scalar quantity nu = 0.5*(lambda + 1/lambda) is commonly
%   used for the real reciprocal pair in planar cases. Here we estimate it
%   from the eigenvalue with largest magnitude.

[V,D] = eig(M);
eigvals = diag(D);

[~, idx] = sort(abs(eigvals), 'descend');
eigvalsSorted = eigvals(idx);

lambdaDom = eigvalsSorted(1);
absDom = abs(lambdaDom);

if abs(lambdaDom) > 0
    nu = 0.5 * (lambdaDom + 1/lambdaDom);
else
    nu = NaN;
end

stab = struct();
stab.eigvals = eigvals;
stab.eigvecs = V;
stab.eigAbsSorted = abs(eigvalsSorted);
stab.lambdaDominant = lambdaDom;
stab.lambdaMin = min(abs(eigvals));
stab.nu = nu;
stab.spectralRadius = max(abs(eigvals));
stab.isApproximatelyStable = abs(stab.spectralRadius - 1) < 1e-6;
end