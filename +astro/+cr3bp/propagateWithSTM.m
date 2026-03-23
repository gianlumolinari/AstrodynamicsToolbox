function out = propagateWithSTM(x0, tspan, mu, opts)
%PROPAGATEWITHSTM Propagate CR3BP state together with the STM.
%
% INPUTS
%   x0    : initial state [6x1]
%   tspan : [t0 tf] or time vector
%   mu    : CR3BP mass parameter
%   opts  : optional odeset options
%
% OUTPUT
%   out : struct with fields
%       .t
%       .Yaug
%       .x
%       .PhiFinal

if nargin < 4 || isempty(opts)
    opts = odeset('RelTol',1e-12, 'AbsTol',1e-12);
end

x0 = x0(:);
Phi0 = eye(6);
Y0 = [x0; Phi0(:)];

[t, Yaug] = ode113(@(t,Y) astro.cr3bp.variationalEOM(t, Y, mu), tspan, Y0, opts);

x = Yaug(:,1:6);
PhiFinal = reshape(Yaug(end,7:end), 6, 6);

out.t = t;
out.Yaug = Yaug;
out.x = x;
out.PhiFinal = PhiFinal;
end