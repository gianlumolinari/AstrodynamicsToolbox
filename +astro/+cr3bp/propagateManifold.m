function out = propagateManifold(xSeed, tManifold, mu)
%PROPAGATEMANIFOLD Propagate one manifold seed in the CR3BP.
%
% INPUTS
%   xSeed     : seed state [6x1]
%   tManifold : final propagation time; use positive for forward,
%               negative for backward propagation
%   mu        : CR3BP mass parameter
%
% OUTPUT
%   out : struct with fields
%       .t
%       .x

xSeed = xSeed(:);

opts.RelTol = 1e-12;
opts.AbsTol = 1e-12;
opts.Solver = 'ode113';

out = astro.propagators.propagate( ...
    @(t,x) astro.cr3bp.eomCR3BP(t, x, mu), ...
    [0 tManifold], xSeed, opts);
end