function n = meanMotion(a, mu)
%MEANMOTION Computes Keplerian mean motion.
%
% INPUTS
%   a  : semi-major axis [km]
%   mu : gravitational parameter [km^3/s^2]
%
% OUTPUT
%   n  : mean motion [rad/s]

if any(a <= 0)
    error('meanMotion is only defined for elliptical orbits with a > 0.');
end

n = sqrt(mu ./ a.^3);
end