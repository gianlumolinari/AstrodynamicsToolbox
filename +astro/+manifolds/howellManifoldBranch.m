function Xseed = howellManifoldBranch(x0, T, mu, nPhase, branchType, sense, epsManifold)
%HOWELLMANIFOLDBRANCH Generate one-sense manifold seeds along a halo orbit.
%
% INPUTS
%   x0          : periodic orbit initial state [6x1]
%   T           : full period
%   mu          : CR3BP mass parameter
%   nPhase      : number of phase samples along orbit
%   branchType  : 'stable' or 'unstable'
%   sense       : +1 or -1
%   epsManifold : perturbation size
%
% OUTPUT
%   Xseed       : [nPhase x 6] manifold seed states
%
% NOTES
%   This is a Ross/Howell-style construction:
%   one branch, one sense, one seed per phase point.

if nargin < 4 || isempty(nPhase)
    nPhase = 80;
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

M = astro.cr3bp.monodromyMatrix(x0, T, mu);
[V,D] = eig(M);
eigvals = diag(D);

switch lower(branchType)
    case 'unstable'
        [~, idx] = max(abs(eigvals));
    case 'stable'
        [~, idx] = min(abs(eigvals));
    otherwise
        error('branchType must be ''stable'' or ''unstable''.');
end

v0 = V(:,idx);

% If the selected branch is slightly complex, use the real part
% for this benchmark-style surface construction.
if norm(imag(v0)) > 1e-8
    warning(['Selected ', branchType, ...
        ' eigenvector is complex; using real part for Howell-style visualization.']);
end
v0 = real(v0);

tGrid = linspace(0, T, nPhase+1);
tGrid(end) = [];

Xseed = zeros(nPhase,6);

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

    nv = norm(vk(1:3));
    if nv < 1e-14
        error('Transported manifold direction norm became too small.');
    end
    vk = vk / nv;

    Xseed(k,:) = (xt + sense*epsManifold*vk).';
end
end