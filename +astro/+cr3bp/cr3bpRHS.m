function dx = cr3bpRHS(~, x, mu)
%CR3BPRHS Right-hand side of the nondimensional CR3BP in synodic frame.
%
% INPUTS
%   x  : state [6x1] = [x; y; z; xd; yd; zd]
%   mu : mass parameter
%
% OUTPUT
%   dx : derivative [6x1]

x = x(:);

xr = x(1);
yr = x(2);
zr = x(3);
xd = x(4);
yd = x(5);
zd = x(6);

r1 = sqrt((xr + mu)^2 + yr^2 + zr^2);
r2 = sqrt((xr - 1 + mu)^2 + yr^2 + zr^2);

xdd = 2*yd + xr ...
    - (1-mu)*(xr + mu)/r1^3 ...
    - mu*(xr - 1 + mu)/r2^3;

ydd = -2*xd + yr ...
    - (1-mu)*yr/r1^3 ...
    - mu*yr/r2^3;

zdd = -(1-mu)*zr/r1^3 ...
    - mu*zr/r2^3;

dx = [xd; yd; zd; xdd; ydd; zdd];
end