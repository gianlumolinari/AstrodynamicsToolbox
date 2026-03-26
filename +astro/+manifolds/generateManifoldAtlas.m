function atlas = generateManifoldAtlas(orbit, mu, type, eps0, timeSpan, opts)
%GENERATEMANIFOLDATLAS Generate a family of manifold branches from one periodic orbit.
%
% INPUTS
%   orbit    : struct with fields .x, .t, .T
%   mu       : CR3BP mass parameter
%   type     : 'stable' or 'unstable'
%   eps0     : seed perturbation magnitude
%   timeSpan : propagation time
%   opts     : optional struct
%
% OUTPUT
%   atlas : struct with fields
%       .type
%       .mu
%       .eps
%       .orbit
%       .eigData
%       .seeds
%       .branches
%       .meta

    if nargin < 6 || isempty(opts)
        opts = struct();
    end

    if ~isfield(opts,'numPhaseSamples'), opts.numPhaseSamples = 40; end
    if ~isfield(opts,'includeBothSigns'), opts.includeBothSigns = true; end
    if ~isfield(opts,'normalizeMode'), opts.normalizeMode = 'position'; end
    if ~isfield(opts,'RelTol'), opts.RelTol = 1e-12; end
    if ~isfield(opts,'AbsTol'), opts.AbsTol = 1e-12; end
    if ~isfield(opts,'Events'), opts.Events = []; end

    eigData = astro.manifolds.computeEigenstructure(orbit, mu, ...
        struct('normalizeMode', opts.normalizeMode, ...
               'RelTol', opts.RelTol, ...
               'AbsTol', opts.AbsTol));

    seeds = astro.manifolds.makeManifoldSeeds(orbit, eigData, type, eps0, ...
        struct('numPhaseSamples', opts.numPhaseSamples, ...
               'includeBothSigns', opts.includeBothSigns));

    nSeeds = numel(seeds);

    branches = struct( ...
        't', {}, ...
        'x', {}, ...
        'type', {}, ...
        'x0', {}, ...
        'C0', {}, ...
        'C', {}, ...
        'CerrMax', {}, ...
        'phaseIndex', {}, ...
        'phaseTime', {}, ...
        'sign', {}, ...
        'xOrbit', {}, ...
        'direction', {});

    for k = 1:nSeeds
        b = astro.manifolds.propagateManifoldBranch( ...
            seeds(k).x0, mu, type, timeSpan, ...
            struct('RelTol', opts.RelTol, ...
                   'AbsTol', opts.AbsTol, ...
                   'Events', opts.Events));

        b.phaseIndex = seeds(k).phaseIndex;
        b.phaseTime  = seeds(k).phaseTime;
        b.sign       = seeds(k).sign;
        b.xOrbit     = seeds(k).xOrbit;
        b.direction  = seeds(k).direction;

        branches(k,1) = b;
    end

    atlas = struct();
    atlas.type = type;
    atlas.mu = mu;
    atlas.eps = eps0;
    atlas.orbit = orbit;
    atlas.eigData = eigData;
    atlas.seeds = seeds;
    atlas.branches = branches;
    atlas.meta.numSeeds = nSeeds;
    atlas.meta.timeSpan = timeSpan;
end