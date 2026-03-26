function seeds = makeManifoldSeeds(orbit, eigData, type, eps0, opts)
%MAKEMANIFOLDSEEDS Generate manifold seed states along a periodic orbit.

    if nargin < 5 || isempty(opts)
        opts = struct();
    end

    if ~isfield(opts,'numPhaseSamples'), opts.numPhaseSamples = 40; end
    if ~isfield(opts,'includeBothSigns'), opts.includeBothSigns = true; end
    if ~isfield(opts,'phaseIndices'), opts.phaseIndices = []; end

    xOrbit = orbit.x;
    tOrbit = orbit.t;
    nOrbit = size(xOrbit,1);

    if isempty(opts.phaseIndices)
        idx = round(linspace(1,nOrbit,opts.numPhaseSamples+1));
        idx(end) = [];
        idx = unique(idx);
    else
        idx = unique(opts.phaseIndices(:).');
    end

    switch lower(type)
        case 'stable'
            dirs = eigData.vStableAlongOrbit;
        case 'unstable'
            dirs = eigData.vUnstableAlongOrbit;
        otherwise
            error('type must be ''stable'' or ''unstable''.');
    end

    if opts.includeBothSigns
        signs = [+1, -1];
    else
        signs = +1;
    end

    seeds = struct('phaseIndex',{},'phaseTime',{},'xOrbit',{}, ...
                   'direction',{},'sign',{},'type',{},'x0',{});

    c = 0;
    for i = 1:numel(idx)
        k = idx(i);
        xk = xOrbit(k,:).';
        vk = dirs(k,:).';

        for s = signs
            c = c + 1;
            seeds(c).phaseIndex = k;
            seeds(c).phaseTime  = tOrbit(k);
            seeds(c).xOrbit     = xk;
            seeds(c).direction  = vk;
            seeds(c).sign       = s;
            seeds(c).type       = type;
            seeds(c).x0         = xk + s * eps0 * vk;
        end
    end
end