function v = visViva(mu, r, a)
%VISVIVA Computes orbital speed from the vis-viva equation.
%
% INPUTS
%   mu : gravitational parameter [km^3/s^2]
%   r  : orbital radius [km]
%   a  : semi-major axis [km]
%
% OUTPUT
%   v  : orbital speed [km/s]

v = sqrt(mu * (2/r - 1/a));
end