function dv = planeChange(v, deltaI)
%PLANECHANGE Computes delta-v for an impulsive plane change.
%
% INPUTS
%   v      : orbital speed [km/s]
%   deltaI : inclination change [rad]
%
% OUTPUT
%   dv     : required delta-v [km/s]

dv = 2 * v * sin(deltaI / 2);
end