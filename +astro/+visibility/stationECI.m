function rECI = stationECI(lat, lon0, h, Re, t, theta0, omegaEarth)
%STATIONECI Ground-station position in ECI using simple Earth rotation.
%
% INPUTS
%   lat        : geodetic latitude [rad]
%   lon0       : longitude at t = 0 [rad]
%   h          : station altitude [km]
%   Re         : Earth radius [km]
%   t          : elapsed time from reference epoch [s]
%   theta0     : Earth rotation angle at t = 0 [rad]
%   omegaEarth : Earth rotation rate [rad/s]
%
% OUTPUT
%   rECI : station position in ECI frame [3x1 km]
%
% MODEL
%   lon(t) = lon0 + theta0 + omegaEarth*t
%
% NOTE
%   This is a simple spherical-Earth, uniform-rotation model.
%   Good for a first access-analysis layer.

lon = lon0 + theta0 + omegaEarth * t;
rECI = astro.visibility.geodeticToECEF(lat, lon, h, Re);
end