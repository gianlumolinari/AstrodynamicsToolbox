function out = flybyFeasibility(vInfIn, vInfOutReq, mu, rpMin, magTol)
%FLYBYFEASIBILITY Checks whether a ballistic flyby can connect two v_infinity vectors.
%
% INPUTS
%   vInfIn     : incoming hyperbolic excess velocity [3x1 km/s]
%   vInfOutReq : required outgoing hyperbolic excess velocity [3x1 km/s]
%   mu         : flyby body gravitational parameter [km^3/s^2]
%   rpMin      : minimum periapsis radius from body center [km]
%   magTol     : allowable magnitude mismatch [km/s] (default: 1e-3)
%
% OUTPUT
%   out : struct with fields
%       .feasible
%       .magIn
%       .magOutReq
%       .magMismatch
%       .magMatch
%       .deltaReqRad
%       .deltaReqDeg
%       .deltaMaxRad
%       .deltaMaxDeg
%       .turnFeasible
%       .rpMin
%       .eccAtRpMin
%
% NOTES
%   A ballistic flyby preserves |v_inf| in the planet-centered frame.
%   Feasibility requires:
%     1) |v_inf^-| ~ |v_inf^+|
%     2) required turning angle <= maximum achievable turning angle

if nargin < 5 || isempty(magTol)
    magTol = 1e-3;
end

vInfIn = vInfIn(:);
vInfOutReq = vInfOutReq(:);

magIn = norm(vInfIn);
magOutReq = norm(vInfOutReq);
magMismatch = abs(magIn - magOutReq);
magMatch = magMismatch <= magTol;

if magIn <= 0 || magOutReq <= 0
    error('Input v_infinity vectors must be nonzero.');
end

cosang = dot(vInfIn, vInfOutReq) / (magIn * magOutReq);
cosang = min(max(cosang, -1), 1);

deltaReqRad = acos(cosang);
deltaReqDeg = rad2deg(deltaReqRad);

% Use incoming magnitude to estimate max turning capability
eccAtRpMin = 1 + rpMin * magIn^2 / mu;
deltaMaxRad = 2 * asin(1 / eccAtRpMin);
deltaMaxDeg = rad2deg(deltaMaxRad);

turnFeasible = deltaReqRad <= deltaMaxRad;
feasible = magMatch && turnFeasible;

out.feasible = feasible;
out.magIn = magIn;
out.magOutReq = magOutReq;
out.magMismatch = magMismatch;
out.magMatch = magMatch;
out.deltaReqRad = deltaReqRad;
out.deltaReqDeg = deltaReqDeg;
out.deltaMaxRad = deltaMaxRad;
out.deltaMaxDeg = deltaMaxDeg;
out.turnFeasible = turnFeasible;
out.rpMin = rpMin;
out.eccAtRpMin = eccAtRpMin;
end