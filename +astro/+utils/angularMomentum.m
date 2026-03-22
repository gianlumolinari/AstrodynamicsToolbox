function h = angularMomentum(r, v)
%ANGULAR MOMENTUM Computes the specific angular momentum vector.
%
% INPUTS
%   r   : position vector [km]
%   v   : velocity vector [km/s]
%
% OUTPUT
%   h   : specific angular momentum vector [km^2/s]

r = r(:);
v = v(:);

h = cross(r, v);
end