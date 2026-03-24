function stab = stabilityIndices(M)
%STABILITYINDICES Compute CR3BP-oriented stability diagnostics from monodromy matrix.
%
% INPUT
%   M : monodromy matrix [6x6]
%
% OUTPUT
%   stab : struct with fields
%       .eigvals
%       .eigvecs
%       .spectralRadius
%       .unstableMultiplier
%       .stableMultiplier
%       .centerPair
%       .unitPair
%       .nu
%       .classification
%       .idxUnstable
%       .idxStable
%       .idxCenter
%       .idxUnit
%       .modulusSorted
%
% NOTES
%   This routine identifies:
%     - unstable multiplier  : largest-modulus eigenvalue
%     - stable multiplier    : smallest-modulus eigenvalue
%     - unit pair            : two eigenvalues closest to +1
%     - center pair          : among the remaining eigenvalues, the pair
%                              whose moduli are closest to 1
%

[V,D] = eig(M);
eigvals = diag(D);
absEig = abs(eigvals);

% Unstable / stable by modulus
[~, idxU] = max(absEig);
[~, idxS] = min(absEig);

lambdaU = eigvals(idxU);
lambdaS = eigvals(idxS);

% Pair closest to +1 in complex plane
distToOne = abs(eigvals - 1);
[~, idxOneSort] = sort(distToOne, 'ascend');
idxUnit = idxOneSort(1:min(2,numel(eigvals)));
unitPair = eigvals(idxUnit);

% Remove unstable/stable candidates, then select remaining pair
allIdx = (1:numel(eigvals)).';
maskRemain = true(size(allIdx));
maskRemain(idxU) = false;
maskRemain(idxS) = false;
remainIdx = allIdx(maskRemain);

centerPair = NaN(2,1);
idxCenter = NaN(2,1);

if numel(remainIdx) >= 2
    remainAbs = abs(eigvals(remainIdx));
    [~, sortRemain] = sort(abs(remainAbs - 1), 'ascend');
    idxCenter = remainIdx(sortRemain(1:2));
    centerPair = eigvals(idxCenter);
end

% Stability index based on unstable multiplier
if abs(lambdaU) > 0
    nu = 0.5 * (lambdaU + 1/lambdaU);
else
    nu = NaN;
end

rho = max(absEig);

if rho < 1 + 1e-8
    cls = 'approximately neutral';
elseif isreal(lambdaU)
    cls = 'real hyperbolic';
else
    cls = 'complex unstable';
end

stab = struct();
stab.eigvals = eigvals;
stab.eigvecs = V;
stab.spectralRadius = rho;
stab.unstableMultiplier = lambdaU;
stab.stableMultiplier = lambdaS;
stab.centerPair = centerPair(:);
stab.unitPair = unitPair(:);
stab.nu = nu;
stab.classification = cls;
stab.idxUnstable = idxU;
stab.idxStable = idxS;
stab.idxCenter = idxCenter(:);
stab.idxUnit = idxUnit(:);
stab.modulusSorted = sort(absEig, 'descend');
end