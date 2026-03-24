function seeds = haloManifoldSeeds(x0, T, mu, nPhase, nCircle, epsManifold)
%HALOMANIFOLDSEEDS Generate halo manifold seeds using 2D stable/unstable subspaces.
%
% INPUTS
%   x0          : periodic halo initial state [6x1]
%   T           : orbit period
%   mu          : CR3BP mass parameter
%   nPhase      : number of phase points along the orbit
%   nCircle     : number of angular samples in each local 2D subspace
%   epsManifold : perturbation radius
%
% OUTPUT
%   seeds : struct array with fields
%       .tPhase
%       .xOrbit
%       .eU1, .eU2
%       .eS1, .eS2
%       .xUnstableRing   [6 x nCircle]
%       .xStableRing     [6 x nCircle]
%
% NOTES
%   The unstable and stable invariant subspaces are constructed from the
%   real and imaginary parts of the complex unstable/stable eigenvectors of
%   the monodromy matrix, then transported along the orbit by the STM.

if nargin < 4 || isempty(nPhase)
    nPhase = 30;
end
if nargin < 5 || isempty(nCircle)
    nCircle = 16;
end
if nargin < 6 || isempty(epsManifold)
    epsManifold = 1e-7;
end

x0 = x0(:);

% ----------------------------------------------------------
% Monodromy eigensystem at initial point
% ----------------------------------------------------------
M = astro.cr3bp.monodromyMatrix(x0, T, mu);
[V,D] = eig(M);
eigvals = diag(D);

[~, idxU] = max(abs(eigvals));
[~, idxS] = min(abs(eigvals));

vU = V(:,idxU);
vS = V(:,idxS);

if abs(imag(vU(1))) < 1e-12 || abs(imag(vS(1))) < 1e-12
    warning('Selected halo manifold eigendirections are nearly real; tube picture may be weak.');
end

eU10 = real(vU);
eU20 = imag(vU);
eS10 = real(vS);
eS20 = imag(vS);

% ----------------------------------------------------------
% Sample orbit phase
% ----------------------------------------------------------
tGrid = linspace(0, T, nPhase+1);
tGrid(end) = [];

theta = linspace(0, 2*pi, nCircle+1);
theta(end) = [];

seeds = repmat(struct( ...
    'tPhase', NaN, ...
    'xOrbit', zeros(6,1), ...
    'eU1', zeros(6,1), ...
    'eU2', zeros(6,1), ...
    'eS1', zeros(6,1), ...
    'eS2', zeros(6,1), ...
    'xUnstableRing', zeros(6,nCircle), ...
    'xStableRing', zeros(6,nCircle)), nPhase, 1);

for k = 1:nPhase
    tk = tGrid(k);

    if abs(tk) < 1e-14
        xt = x0;
        Phi_t0 = eye(6);
    else
        prop = astro.cr3bp.propagateWithSTM( ...
            x0, [0 tk], mu, odeset('RelTol',1e-12,'AbsTol',1e-12));
        xt = prop.x(end,:).';
        Phi_t0 = prop.PhiFinal;
    end

    % Transport subspace bases
    eU1 = Phi_t0 * eU10;
    eU2 = Phi_t0 * eU20;
    eS1 = Phi_t0 * eS10;
    eS2 = Phi_t0 * eS20;

    % Orthonormalize approximately using 3D position part
    [eU1, eU2] = localOrthonormalize(eU1, eU2);
    [eS1, eS2] = localOrthonormalize(eS1, eS2);

    XU = zeros(6,nCircle);
    XS = zeros(6,nCircle);

    for j = 1:nCircle
        du = cos(theta(j))*eU1 + sin(theta(j))*eU2;
        ds = cos(theta(j))*eS1 + sin(theta(j))*eS2;

        XU(:,j) = xt + epsManifold * du;
        XS(:,j) = xt + epsManifold * ds;
    end

    seeds(k).tPhase = tk;
    seeds(k).xOrbit = xt;
    seeds(k).eU1 = eU1;
    seeds(k).eU2 = eU2;
    seeds(k).eS1 = eS1;
    seeds(k).eS2 = eS2;
    seeds(k).xUnstableRing = XU;
    seeds(k).xStableRing = XS;
end

end

function [e1, e2] = localOrthonormalize(v1, v2)
% Use position-part Gram-Schmidt for cleaner geometric tube directions.

e1 = v1;
n1 = norm(e1(1:3));
if n1 < 1e-14
    error('Local manifold basis vector norm became too small.');
end
e1 = e1 / n1;

alpha = dot(v2(1:3), e1(1:3));
e2 = v2 - alpha * e1;

n2 = norm(e2(1:3));
if n2 < 1e-14
    % fallback if nearly collinear
    e2 = v2;
    n2 = norm(e2(1:3));
    if n2 < 1e-14
        error('Second local manifold basis vector norm became too small.');
    end
end
e2 = e2 / n2;
end