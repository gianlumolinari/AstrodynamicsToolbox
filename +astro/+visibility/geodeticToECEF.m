function rECEF = geodeticToECEF(lat, lon, h, Re)
%GEODETICTOECEF Simple spherical Earth geodetic to ECEF conversion.
%
% INPUTS
%   lat : geodetic latitude [rad]
%   lon : longitude [rad]
%   h   : altitude above spherical Earth [km]
%   Re  : Earth radius [km]
%
% OUTPUT
%   rECEF : position in ECEF frame [3x1 km]
%
% NOTE
%   This first version assumes a spherical Earth.

rho = Re + h;

x = rho * cos(lat) * cos(lon);
y = rho * cos(lat) * sin(lon);
z = rho * sin(lat);

rECEF = [x; y; z];
end