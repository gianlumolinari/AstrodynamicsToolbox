function eps = specificEnergy(r, v, mu)
%SPECIFICENERGY Computes the specific orbital energy.
%
% INPUTS
%   r   : position vector [km]
%   v   : velocity vector [km/s]
%   mu  : gravitational parameter [km^3/s^2]
%
% OUTPUT
%   eps : specific orbital energy [km^2/s^2]

r = r(:);
v = v(:);

eps = 0.5 * dot(v, v) - mu / norm(r);
end