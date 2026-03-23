function T = orbitalPeriod(a, mu)
%ORBITALPERIOD Orbital period for a Keplerian ellipse.
%
% INPUTS
%   a  : semi-major axis [km]
%   mu : gravitational parameter [km^3/s^2]
%
% OUTPUT
%   T  : orbital period [s]
%
% NOTE
%   Valid for elliptical orbits only (a > 0).

if any(a <= 0)
    error('orbitalPeriod is only defined for elliptical orbits with a > 0.');
end

T = 2*pi*sqrt(a.^3 ./ mu);
end