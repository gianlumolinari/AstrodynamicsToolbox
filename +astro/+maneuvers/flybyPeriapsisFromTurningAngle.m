function out = flybyPeriapsisFromTurningAngle(muPlanet, vInf, delta)
%FLYBYPERIAPSISFROMTURNINGANGLE Computes periapsis radius needed for a flyby turn.
%
% INPUTS
%   muPlanet : planet gravitational parameter [km^3/s^2]
%   vInf     : hyperbolic excess speed or vector [km/s]
%   delta    : turning angle [rad]
%
% OUTPUT
%   out : struct with fields
%       .vInf
%       .deltaRad
%       .deltaDeg
%       .e
%       .rp
%
% Relations:
%   delta = 2*asin(1/e)   =>   e = 1/sin(delta/2)
%   rp = mu*(e - 1)/vInf^2

if isvector(vInf) && numel(vInf) == 3
    vInfMag = norm(vInf);
else
    vInfMag = vInf;
end

e = 1 / sin(delta / 2);
rp = muPlanet * (e - 1) / vInfMag^2;

out.vInf = vInfMag;
out.deltaRad = delta;
out.deltaDeg = rad2deg(delta);
out.e = e;
out.rp = rp;
end