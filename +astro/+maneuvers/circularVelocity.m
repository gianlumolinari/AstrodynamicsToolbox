function v = circularVelocity(mu, r)
%CIRCULARVELOCITY Circular orbital speed.
%
% INPUTS
%   mu : gravitational parameter [km^3/s^2]
%   r  : orbital radius [km]
%
% OUTPUT
%   v  : circular speed [km/s]

v = sqrt(mu ./ r);
end