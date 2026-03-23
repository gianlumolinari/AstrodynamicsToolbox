function dx = eomCowell(~, x, body, pert)
%EOMCOWELL Cowell equations with optional perturbations.
%
% INPUTS
%   x    : state [6x1] = [r; v]
%   body : struct with fields mu, radius, J2
%   pert : struct with optional fields:
%          .useJ2
%          .useDrag
%          .Cd
%          .AoverM
%          .rho0
%          .H
%          .omegaBody
%
% OUTPUT
%   dx : state derivative [6x1]

r = x(1:3);
v = x(4:6);

rmag = norm(r);
a = -body.mu * r / rmag^3;

if nargin >= 4 && ~isempty(pert)
    if isfield(pert, 'useJ2') && pert.useJ2
        a = a + astro.perturbations.accelJ2(r, body.mu, body.J2, body.radius);
    end

    if isfield(pert, 'useDrag') && pert.useDrag
        a = a + astro.perturbations.accelDrag( ...
            r, v, pert.Cd, pert.AoverM, pert.rho0, pert.H, body.radius, pert.omegaBody);
    end
end

dx = [v; a];
end