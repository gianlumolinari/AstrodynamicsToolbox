function out = propagate(dynFun, tspan, x0, opts)
%PROPAGATE Generic numerical propagator wrapper.
%
% INPUTS
%   dynFun : function handle of the form @(t,x) ...
%   tspan  : integration interval [t0 tf]
%   x0     : initial state vector
%   opts   : optional struct with fields:
%       .RelTol
%       .AbsTol
%       .Solver   ('ode45' or 'ode113')
%
% OUTPUT
%   out : struct with fields
%       .t
%       .x

if nargin < 4
    opts = struct();
end

if ~isfield(opts, 'RelTol')
    opts.RelTol = 1e-12;
end

if ~isfield(opts, 'AbsTol')
    opts.AbsTol = 1e-12;
end

if ~isfield(opts, 'Solver')
    opts.Solver = 'ode113';
end

odeOpts = odeset('RelTol', opts.RelTol, 'AbsTol', opts.AbsTol);

switch lower(opts.Solver)
    case 'ode45'
        [t, x] = ode45(dynFun, tspan, x0, odeOpts);
    case 'ode113'
        [t, x] = ode113(dynFun, tspan, x0, odeOpts);
    otherwise
        error('Unsupported solver. Use ''ode45'' or ''ode113''.');
end

out.t = t;
out.x = x;
end