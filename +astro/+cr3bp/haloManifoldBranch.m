function Xseed = haloManifoldBranch(x0, T, mu, nPhase, branchType, sense, epsManifold)
%HALOMANIFOLDBRANCH Generate Ross-style halo manifold seeds.
%
% INPUTS
%   x0          : periodic orbit initial state [6x1]
%   T           : orbit period
%   mu          : CR3BP mass parameter
%   nPhase      : number of phase points along the orbit
%   branchType  : 'stable' or 'unstable'
%   sense       : +1 or -1
%   epsManifold : perturbation size
%
% OUTPUT
%   Xseed       : [nPhase x 6] manifold initial conditions
%
% NOTES
%   This uses one eigendirection and one sign, which is the right setup
%   for Ross-style manifold surface plots.

if nargin < 4 || isempty(nPhase)
    nPhase = 60;
end
if nargin < 5 || isempty(branchType)
    branchType = 'stable';
end
if nargin < 6 || isempty(sense)
    sense = 1;
end
if nargin < 7 || isempty(epsManifold)
    epsManifold = 1e-6;
end

x0 = x0(:);

% Monodromy eigensystem
M = astro.cr3bp.monodromyMatrix(x0, T, mu);
[V,D] = eig(M);
eigvals = diag(D);

% Choose branch
switch lower(branchType)
    case 'unstable'
        [~, idx] = max(abs(eigvals));
    case 'stable'
        [~, idx] = min(abs(eigvals));
    otherwise
        error('branchType must be ''stable'' or ''unstable''.');
end

v0 = V(:,idx);

% This Ross-style construction expects a real branch
if norm(imag(v0)) > 1e-8
    error(['Selected ', branchType, ...
        ' eigenvector is complex. Choose a different orbit/member, or use a different manifold construction.']);
end

v0 = real(v0);

% Phase sampling
tGrid = linspace(0, T, nPhase+1);
tGrid(end) = [];

Xseed = zeros(nPhase, 6);

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

    vk = Phi_t0 * v0;

    % normalize by position part
    nv = norm(vk(1:3));
    if nv < 1e-14
        error('Transported manifold direction norm became too small.');
    end
    vk = vk / nv;

    Xseed(k,:) = (xt + sense*epsManifold*vk).';
end
end