function out = azimuthElevationRange(rSc, rStation, lat, lon)
%AZIMUTHELEVATIONRANGE Computes topocentric azimuth, elevation, and range.
%
% INPUTS
%   rSc      : spacecraft position [3x1 km]
%   rStation : station position [3x1 km]
%   lat      : station latitude [rad]
%   lon      : station longitude [rad]
%
% OUTPUT
%   out : struct with fields
%       .azimuthRad
%       .azimuthDeg
%       .elevationRad
%       .elevationDeg
%       .range
%       .sez

sez = astro.visibility.topocentricSEZ(rSc, rStation, lat, lon);

S = sez(1);
E = sez(2);
Z = sez(3);

range = norm(sez);
elevationRad = asin(Z / range);

azimuthRad = atan2(E, -S);
if azimuthRad < 0
    azimuthRad = azimuthRad + 2*pi;
end

out.azimuthRad = azimuthRad;
out.azimuthDeg = rad2deg(azimuthRad);
out.elevationRad = elevationRad;
out.elevationDeg = rad2deg(elevationRad);
out.range = range;
out.sez = sez;
end