function rSOI = sphereOfInfluence(aPlanet, muPlanet, muSun)
%SPHEREOFINFLUENCE Computes Laplace sphere of influence radius.
%
% INPUTS
%   aPlanet  : heliocentric semimajor axis of planet [km]
%   muPlanet : planet gravitational parameter [km^3/s^2]
%   muSun    : Sun gravitational parameter [km^3/s^2]
%
% OUTPUT
%   rSOI     : sphere of influence radius [km]

rSOI = aPlanet * (muPlanet / muSun)^(2/5);
end