function C = jacobiConstant(x, mu)
%JACOBICONSTANT Computes the Jacobi constant in the CR3BP.
%
% INPUTS
%   x  : state [6x1] or [Nx6]
%   mu : mass parameter
%
% OUTPUT
%   C  : Jacobi constant

if size(x,2) ~= 6
    x = x.';
end

xr = x(:,1);
yr = x(:,2);
zr = x(:,3);
xd = x(:,4);
yd = x(:,5);
zd = x(:,6);

r1 = sqrt((xr + mu).^2 + yr.^2 + zr.^2);
r2 = sqrt((xr - 1 + mu).^2 + yr.^2 + zr.^2);

Omega = 0.5*(xr.^2 + yr.^2) + (1-mu)./r1 + mu./r2;
v2 = xd.^2 + yd.^2 + zd.^2;

C = 2*Omega - v2;
end