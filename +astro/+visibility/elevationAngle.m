function out = elevationAngle(rSc, rStation)
%ELEVATIONANGLE Computes spacecraft elevation angle from a ground station.
%
% INPUTS
%   rSc      : spacecraft position in Earth-centered inertial frame [3x1 km]
%   rStation : station position in same frame [3x1 km]
%
% OUTPUT
%   out : struct with fields
%       .elevationRad
%       .elevationDeg
%       .range
%       .losHat
%       .zenithHat
%
% NOTE
%   Uses a simple spherical Earth geometric elevation definition.

rSc = rSc(:);
rStation = rStation(:);

rho = rSc - rStation;
range = norm(rho);

if range <= 0
    error('Spacecraft and station positions must be distinct.');
end

losHat = rho / range;
zenithHat = rStation / norm(rStation);

elevationRad = asin(dot(losHat, zenithHat));
elevationDeg = rad2deg(elevationRad);

out.elevationRad = elevationRad;
out.elevationDeg = elevationDeg;
out.range = range;
out.losHat = losHat;
out.zenithHat = zenithHat;
end