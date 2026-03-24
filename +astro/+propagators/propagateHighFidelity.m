function out = propagateHighFidelity(x0, et0, tf, config)
%PROPAGATEHIGHFIDELITY Propagate heliocentric trajectory with high-fidelity model.
%
% INPUTS
%   x0     : initial state wrt Sun [km; km/s]
%   et0    : initial SPICE ET
%   tf     : final propagation time since et0 [s]
%   config : configuration struct
%
% OUTPUT
%   out.t  : propagation times [s]
%   out.x  : propagated states

    if ~isfield(config,'RelTol'); config.RelTol = 1e-12; end
    if ~isfield(config,'AbsTol'); config.AbsTol = 1e-12; end

    opts = odeset('RelTol', config.RelTol, 'AbsTol', config.AbsTol);

    rhs = @(t,x) astro.propagators.eomNBodySpiceSRP(t, x, config);

    [t, x] = ode113(rhs, [0 tf], x0, opts);

    out.t = t;
    out.x = x;
end