function out = hyperbolicPeriapsisRadius(muPlanet, vInf, rp)
%HYPERBOLICPERIAPSISRADIUS Computes hyperbolic parameters from periapsis radius.
%
% INPUTS
%   muPlanet : planet gravitational parameter [km^3/s^2]
%   vInf     : hyperbolic excess speed or vector [km/s]
%   rp       : periapsis radius from planet center [km]
%
% OUTPUT
%   out : struct with fields
%       .vInf
%       .rp
%       .e
%       .a
%       .vp
%
% Relations:
%   e = 1 + rp*vInf^2/mu
%   a = -mu/vInf^2
%   vp = sqrt(vInf^2 + 2*mu/rp)

if isvector(vInf) && numel(vInf) == 3
    vInfMag = norm(vInf);
else
    vInfMag = vInf;
end

e = 1 + rp * vInfMag^2 / muPlanet;
a = -muPlanet / vInfMag^2;
vp = sqrt(vInfMag^2 + 2*muPlanet / rp);

out.vInf = vInfMag;
out.rp = rp;
out.e = e;
out.a = a;
out.vp = vp;
end