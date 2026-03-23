function Omega = effectivePotential(x, y, z, mu)
%EFFECTIVEPOTENTIAL Effective potential Omega in the CR3BP.
%
% INPUTS
%   x, y, z : coordinates (scalars or arrays)
%   mu      : mass parameter
%
% OUTPUT
%   Omega   : effective potential

r1 = sqrt((x + mu).^2 + y.^2 + z.^2);
r2 = sqrt((x - 1 + mu).^2 + y.^2 + z.^2);

Omega = 0.5*(x.^2 + y.^2) + (1-mu)./r1 + mu./r2;
end