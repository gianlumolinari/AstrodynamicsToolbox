function seeds = manifoldSeeds(x0, T, mu, nPhase, epsManifold)
%MANIFOLDSEEDS Generate stable/unstable manifold seeds along a periodic orbit.
%
% INPUTS
%   x0          : periodic-orbit initial state [6x1]
%   T           : orbit period
%   mu          : CR3BP mass parameter
%   nPhase      : number of phase points along the orbit
%   epsManifold : perturbation magnitude
%
% OUTPUT
%   seeds : struct array with fields
%       .tPhase
%       .xOrbit
%       .vStableLocal
%       .vUnstableLocal
%       .xStablePlus
%       .xStableMinus
%       .xUnstablePlus
%       .xUnstableMinus

if nargin < 4 || isempty(nPhase)
    nPhase = 20;
end
if nargin < 5 || isempty(epsManifold)
    epsManifold = 1e-6;
end

x0 = x0(:);

% Monodromy and eigenstructure at initial point
M = astro.cr3bp.monodromyMatrix(x0, T, mu);
[V,D] = eig(M);
eigvals = diag(D);

[~, idxU] = max(abs(eigvals));
[~, idxS] = min(abs(eigvals));

vU0 = real(V(:,idxU));
vS0 = real(V(:,idxS));

% Sample one full orbit, excluding duplicate endpoint
tGrid = linspace(0, T, nPhase+1);
tGrid(end) = [];

seeds = repmat(struct( ...
    'tPhase', NaN, ...
    'xOrbit', zeros(6,1), ...
    'vStableLocal', zeros(6,1), ...
    'vUnstableLocal', zeros(6,1), ...
    'xStablePlus', zeros(6,1), ...
    'xStableMinus', zeros(6,1), ...
    'xUnstablePlus', zeros(6,1), ...
    'xUnstableMinus', zeros(6,1)), nPhase, 1);

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

    vU = Phi_t0 * vU0;
    vS = Phi_t0 * vS0;

    % Normalize using position-part norm
    nU = norm(vU(1:3));
    nS = norm(vS(1:3));

    if nU < 1e-14 || nS < 1e-14
        error('Local manifold eigenvector norm became too small.');
    end

    vU = vU / nU;
    vS = vS / nS;

    seeds(k).tPhase = tk;
    seeds(k).xOrbit = xt;
    seeds(k).vStableLocal = vS;
    seeds(k).vUnstableLocal = vU;

    seeds(k).xStablePlus = xt + epsManifold * vS;
    seeds(k).xStableMinus = xt - epsManifold * vS;
    seeds(k).xUnstablePlus = xt + epsManifold * vU;
    seeds(k).xUnstableMinus = xt - epsManifold * vU;
end

end