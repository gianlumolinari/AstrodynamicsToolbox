function v = escapeVelocity(mu, r)
%ESCAPEVELOCITY Escape speed at radius r.
%
% INPUTS
%   mu : gravitational parameter [km^3/s^2]
%   r  : orbital radius [km]
%
% OUTPUT
%   v  : escape speed [km/s]

v = sqrt(2 * mu ./ r);
end