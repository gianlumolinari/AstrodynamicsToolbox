function out = hyperbolicTurningAngle(muPlanet, rp, vInf)
%HYPERBOLICTURNINGANGLE Computes hyperbolic flyby turning angle.
%
% INPUTS
%   muPlanet : planet gravitational parameter [km^3/s^2]
%   rp       : hyperbolic periapsis radius from planet center [km]
%   vInf     : hyperbolic excess speed or vector [km/s]
%
% OUTPUT
%   out : struct with fields
%       .vInf
%       .rp
%       .e
%       .deltaRad
%       .deltaDeg
%
% Relations:
%   e = 1 + rp*vInf^2/mu
%   delta = 2*asin(1/e)

if isvector(vInf) && numel(vInf) == 3
    vInfMag = norm(vInf);
else
    vInfMag = vInf;
end

e = 1 + rp * vInfMag^2 / muPlanet;
deltaRad = 2 * asin(1 / e);
deltaDeg = rad2deg(deltaRad);

out.vInf = vInfMag;
out.rp = rp;
out.e = e;
out.deltaRad = deltaRad;
out.deltaDeg = deltaDeg;
end