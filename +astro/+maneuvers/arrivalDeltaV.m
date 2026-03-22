function out = arrivalDeltaV(muPlanet, rParking, vInf)
%ARRIVALDELTAV Computes impulsive capture delta-v into circular orbit.
%
% INPUTS
%   muPlanet : planet gravitational parameter [km^3/s^2]
%   rParking : radius of circular capture orbit from planet center [km]
%   vInf     : hyperbolic excess speed or vector [km/s]
%
% OUTPUT
%   out : struct with fields
%       .vInf
%       .vCirc
%       .vPeriHyp
%       .deltaV

if isvector(vInf) && numel(vInf) == 3
    vInfMag = norm(vInf);
else
    vInfMag = vInf;
end

vCirc = sqrt(muPlanet / rParking);
vPeriHyp = sqrt(vInfMag^2 + 2*muPlanet / rParking);
deltaV = vPeriHyp - vCirc;

out.vInf = vInfMag;
out.vCirc = vCirc;
out.vPeriHyp = vPeriHyp;
out.deltaV = deltaV;
end