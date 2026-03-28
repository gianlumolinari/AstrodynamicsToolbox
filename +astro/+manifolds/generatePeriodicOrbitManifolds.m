function man = generatePeriodicOrbitManifolds(orbit, mu, opts)
%GENERATEPERIODICORBITMANIFOLDS Unified manifold generation for any periodic orbit.
%
% INPUTS
%   orbit : generic periodic-orbit struct
%   mu    : CR3BP mass parameter
%   opts  : struct with fields
%       .eps0
%       .tfStable
%       .tfUnstable
%       .numPhaseSamples
%       .includeBothSigns
%       .normalizeMode
%       .RelTol
%       .AbsTol
%       .validateOrbit
%       .verbose
%
% OUTPUT
%   man : struct with fields
%       .orbit
%       .validation
%       .stable
%       .unstable
%       .maxJacobiDriftStable
%       .maxJacobiDriftUnstable

if nargin < 3 || isempty(opts)
    opts = struct();
end

if ~isfield(opts,'eps0'), opts.eps0 = 2e-5; end
if ~isfield(opts,'tfStable'), opts.tfStable = 5.0; end
if ~isfield(opts,'tfUnstable'), opts.tfUnstable = 5.0; end
if ~isfield(opts,'numPhaseSamples'), opts.numPhaseSamples = 25; end
if ~isfield(opts,'includeBothSigns'), opts.includeBothSigns = true; end
if ~isfield(opts,'normalizeMode'), opts.normalizeMode = 'position'; end
if ~isfield(opts,'RelTol'), opts.RelTol = 1e-12; end
if ~isfield(opts,'AbsTol'), opts.AbsTol = 1e-12; end
if ~isfield(opts,'validateOrbit'), opts.validateOrbit = true; end
if ~isfield(opts,'verbose'), opts.verbose = true; end

validation = [];
if opts.validateOrbit
    validation = astro.manifolds.validatePeriodicOrbitForManifolds( ...
        orbit, mu, struct('verbose', opts.verbose));
    if ~validation.isAcceptable
        error('Periodic orbit failed manifold validation. Fix orbit before generating manifolds.');
    end
end

atlasU = astro.manifolds.generateManifoldAtlas(orbit, mu, 'unstable', opts.eps0, opts.tfUnstable, ...
    struct('numPhaseSamples', opts.numPhaseSamples, ...
           'includeBothSigns', opts.includeBothSigns, ...
           'normalizeMode', opts.normalizeMode, ...
           'RelTol', opts.RelTol, ...
           'AbsTol', opts.AbsTol));

atlasS = astro.manifolds.generateManifoldAtlas(orbit, mu, 'stable', opts.eps0, opts.tfStable, ...
    struct('numPhaseSamples', opts.numPhaseSamples, ...
           'includeBothSigns', opts.includeBothSigns, ...
           'normalizeMode', opts.normalizeMode, ...
           'RelTol', opts.RelTol, ...
           'AbsTol', opts.AbsTol));

atlasU = astro.manifolds.sortTubeStrands(atlasU, 'signphase');
atlasS = astro.manifolds.sortTubeStrands(atlasS, 'signphase');

man = struct();
man.orbit = orbit;
man.validation = validation;
man.unstable = atlasU;
man.stable = atlasS;
man.maxJacobiDriftUnstable = max(arrayfun(@(b) b.CerrMax, atlasU.branches));
man.maxJacobiDriftStable   = max(arrayfun(@(b) b.CerrMax, atlasS.branches));

if opts.verbose
    fprintf('\nUnified manifold summary\n');
    fprintf('------------------------\n');
    fprintf('Unstable branches: %d\n', numel(atlasU.branches));
    fprintf('Stable branches  : %d\n', numel(atlasS.branches));
    fprintf('Max Jacobi drift (unstable): %.3e\n', man.maxJacobiDriftUnstable);
    fprintf('Max Jacobi drift (stable)  : %.3e\n', man.maxJacobiDriftStable);
end
end