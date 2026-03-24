function dY = variationalEOM(~, Y, mu)
%VARIATIONALEOM CR3BP equations plus STM dynamics.
%
% State packing:
%   Y = [x;y;z;xd;yd;zd; reshape(Phi,36,1)]
%
% Output:
%   dY = [state_dot; reshape(Phi_dot,36,1)]

x  = Y(1);
y  = Y(2);
z  = Y(3);
xd = Y(4);
yd = Y(5);
zd = Y(6);

r1 = sqrt((x + mu)^2 + y^2 + z^2);
r2 = sqrt((x - 1 + mu)^2 + y^2 + z^2);

% State equations
ddx = 2*yd + x - (1-mu)*(x+mu)/r1^3 - mu*(x-1+mu)/r2^3;
ddy = -2*xd + y - (1-mu)*y/r1^3 - mu*y/r2^3;
ddz = -(1-mu)*z/r1^3 - mu*z/r2^3;

% Second derivatives of potential
Uxx = 1 ...
    - (1-mu)/r1^3 - mu/r2^3 ...
    + 3*(1-mu)*(x+mu)^2/r1^5 ...
    + 3*mu*(x-1+mu)^2/r2^5;

Uyy = 1 ...
    - (1-mu)/r1^3 - mu/r2^3 ...
    + 3*(1-mu)*y^2/r1^5 ...
    + 3*mu*y^2/r2^5;

Uzz = -(1-mu)/r1^3 - mu/r2^3 ...
    + 3*(1-mu)*z^2/r1^5 ...
    + 3*mu*z^2/r2^5;

Uxy = 3*(1-mu)*(x+mu)*y/r1^5 + 3*mu*(x-1+mu)*y/r2^5;
Uxz = 3*(1-mu)*(x+mu)*z/r1^5 + 3*mu*(x-1+mu)*z/r2^5;
Uyz = 3*(1-mu)*y*z/r1^5 + 3*mu*y*z/r2^5;

A = [0   0   0   1  0  0;
     0   0   0   0  1  0;
     0   0   0   0  0  1;
     Uxx Uxy Uxz 0  2  0;
     Uxy Uyy Uyz -2 0  0;
     Uxz Uyz Uzz 0  0  0];

Phi = reshape(Y(7:end), 6, 6);
PhiDot = A * Phi;

dY = [xd; yd; zd; ddx; ddy; ddz; PhiDot(:)];
end