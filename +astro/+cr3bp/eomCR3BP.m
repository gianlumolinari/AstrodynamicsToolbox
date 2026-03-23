function dx = eomCR3BP(~, x, mu)
%EOMCR3BP Equations of motion for the nondimensional circular restricted
% three-body problem in the synodic frame.
%
% INPUTS
%   x  : state [6x1] = [x; y; z; xd; yd; zd]
%   mu : mass parameter, e.g. Earth-Moon mu ~ 0.01215
%
% OUTPUT
%   dx : state derivative [6x1]

x = x(:);

xr  = x(1);
yr  = x(2);
zr  = x(3);
xd  = x(4);
yd  = x(5);
zd  = x(6);

r1 = sqrt((xr + mu)^2 + yr^2 + zr^2);
r2 = sqrt((xr - 1 + mu)^2 + yr^2 + zr^2);

ddx = 2*yd + xr ...
    - (1-mu)*(xr + mu)/r1^3 ...
    - mu*(xr - 1 + mu)/r2^3;

ddy = -2*xd + yr ...
    - (1-mu)*yr/r1^3 ...
    - mu*yr/r2^3;

ddz = -(1-mu)*zr/r1^3 ...
      - mu*zr/r2^3;

dx = [xd; yd; zd; ddx; ddy; ddz];
end