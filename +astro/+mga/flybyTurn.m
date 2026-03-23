function out = flybyTurn(vInfIn, mu, rp, turnSign)
%FLYBYTURN Planar gravity-assist turn of a hyperbolic excess velocity vector.
%
% INPUTS
%   vInfIn   : incoming hyperbolic excess velocity [3x1 km/s]
%   mu       : flyby body gravitational parameter [km^3/s^2]
%   rp       : periapsis radius from body center [km]
%   turnSign : +1 or -1 for left/right turn sense
%
% OUTPUT
%   out : struct with fields
%       .vInfIn
%       .vInfOut
%       .vInfMag
%       .deltaRad
%       .deltaDeg
%       .ecc
%       .rp
%
% NOTES
%   Planar model: the outgoing vector is produced by rotating the incoming
%   vector in the xy-plane by the turning angle.

if nargin < 4 || isempty(turnSign)
    turnSign = +1;
end

vInfIn = vInfIn(:);
vInfMag = norm(vInfIn);

if vInfMag <= 0
    error('vInfIn must have nonzero magnitude.');
end
if rp <= 0
    error('rp must be positive.');
end
if abs(turnSign) ~= 1
    error('turnSign must be +1 or -1.');
end

ecc = 1 + rp * vInfMag^2 / mu;
deltaRad = 2 * asin(1 / ecc);
deltaRad = turnSign * deltaRad;
deltaDeg = rad2deg(deltaRad);

R = [cos(deltaRad), -sin(deltaRad), 0;
     sin(deltaRad),  cos(deltaRad), 0;
     0,              0,             1];

vInfOut = R * vInfIn;

out.vInfIn = vInfIn;
out.vInfOut = vInfOut;
out.vInfMag = vInfMag;
out.deltaRad = deltaRad;
out.deltaDeg = deltaDeg;
out.ecc = ecc;
out.rp = rp;
end