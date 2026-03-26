function branch = propagateManifoldBranch(x0, mu, type, timeSpan, opts)
%PROPAGATEMANIFOLDBRANCH Propagate one stable or unstable manifold branch.

    if nargin < 5 || isempty(opts)
        opts = struct();
    end

    if ~isfield(opts,'RelTol'), opts.RelTol = 1e-12; end
    if ~isfield(opts,'AbsTol'), opts.AbsTol = 1e-12; end
    if ~isfield(opts,'Events'), opts.Events = []; end

    if isscalar(timeSpan)
        switch lower(type)
            case 'unstable'
                tspan = [0, abs(timeSpan)];
            case 'stable'
                tspan = [0, -abs(timeSpan)];
            otherwise
                error('type must be ''stable'' or ''unstable''.');
        end
    else
        tspan = timeSpan(:).';
    end

    odeOpts = odeset('RelTol',opts.RelTol,'AbsTol',opts.AbsTol);
    if ~isempty(opts.Events)
        odeOpts = odeset(odeOpts,'Events',opts.Events);
    end

    rhs = @(t,x) astro.cr3bp.eomCR3BP(t, x, mu);
    [t, x] = ode113(rhs, tspan, x0, odeOpts);

    C = zeros(size(t));
    for k = 1:numel(t)
        C(k) = astro.cr3bp.jacobiConstant(x(k,:).', mu);
    end

    C0 = C(1);

    branch = struct();
    branch.t = t;
    branch.x = x;
    branch.type = type;
    branch.x0 = x0;
    branch.C0 = C0;
    branch.C = C;
    branch.CerrMax = max(abs(C - C0));
end