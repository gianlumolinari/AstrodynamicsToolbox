function sez = topocentricSEZ(rSc, rStation, lat, lon)
%TOPOCENTRICSEZ Converts spacecraft LOS vector to local SEZ frame.
%
% INPUTS
%   rSc      : spacecraft position in ECI-like frame [3x1 km]
%   rStation : station position in same frame [3x1 km]
%   lat      : station latitude [rad]
%   lon      : station longitude [rad]
%
% OUTPUT
%   sez : LOS vector in SEZ frame [3x1 km]
%
% NOTE
%   This uses a simple local South-East-Zenith frame.

rho = rSc(:) - rStation(:);

slat = sin(lat); clat = cos(lat);
slon = sin(lon); clon = cos(lon);

R = [ slat*clon,  slat*slon, -clat;
     -slon,       clon,       0;
      clat*clon,  clat*slon,  slat ];

sez = R * rho;
end