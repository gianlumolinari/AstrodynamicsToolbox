function dx = eomTwoBody(~, x, mu)
%EOMTWOBODY Two-body equations of motion.
%
% INPUTS
%   x   : state vector [rx; ry; rz; vx; vy; vz]
%   mu  : gravitational parameter [km^3/s^2]
%
% OUTPUT
%   dx  : time derivative of the state vector

r = x(1:3);
v = x(4:6);

rnorm = norm(r);
a = -mu / rnorm^3 * r;

dx = [v; a];
end