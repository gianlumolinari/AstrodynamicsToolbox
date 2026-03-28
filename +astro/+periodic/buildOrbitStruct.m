function orbit = buildOrbitStruct(state0, T, mu, meta, opts)
%BUILDORBITSTRUCT Build a generic periodic-orbit struct for manifold work.
%
% INPUTS
%   state0 : 6x1 periodic-orbit initial state
%   T      : period
%   mu     : CR3BP mass parameter
%   meta   : optional struct with fields:
%            .family, .libr, .system, .source, .dimension
%   opts   : optional integrator settings
%
% OUTPUT
%   orbit  : struct with standard fields for manifold generation

if nargin < 4 || isempty(meta)
    meta = struct();
end
if nargin < 5 || isempty(opts)
    opts = struct();
end

if ~isfield(opts,'RelTol'), opts.RelTol = 1e-12; end
if ~isfield(opts,'AbsTol'), opts.AbsTol = 1e-12; end
if ~isfield(opts,'verbose'), opts.verbose = true; end

prop = astro.cr3bp.propagateWithSTM( ...
    state0, [0 T], mu, odeset('RelTol',opts.RelTol,'AbsTol',opts.AbsTol));

orbit = struct();
orbit.state0 = state0;
orbit.T = T;
orbit.mu = mu;
orbit.t = prop.t(:);
orbit.x = prop.x(:,1:6);
orbit.monodromy = prop.PhiFinal;
orbit.closureError = norm(orbit.x(end,:).' - state0);

if ~isfield(meta,'family'), meta.family = ''; end
if ~isfield(meta,'libr'), meta.libr = []; end
if ~isfield(meta,'system'), meta.system = ''; end
if ~isfield(meta,'source'), meta.source = ''; end
if ~isfield(meta,'dimension'), meta.dimension = '3D'; end

orbit.family = meta.family;
orbit.libr = meta.libr;
orbit.system = meta.system;
orbit.source = meta.source;
orbit.dimension = meta.dimension;

if opts.verbose
    fprintf('Orbit closure error after one period: %.3e\n', orbit.closureError);
end
end