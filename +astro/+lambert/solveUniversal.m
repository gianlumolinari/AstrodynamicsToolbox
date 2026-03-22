function sol = solveUniversal(r1, r2, tof, mu, longWay)
%SOLVEUNIVERSAL Solves the 0-rev Lambert problem using universal variables.
%
% INPUTS
%   r1      : initial position vector [km]
%   r2      : final position vector [km]
%   tof     : time of flight [s]
%   mu      : gravitational parameter [km^3/s^2]
%   longWay : logical flag
%             false -> short-way transfer
%             true  -> long-way transfer
%
% OUTPUT
%   sol : struct with fields
%       .v1
%       .v2
%       .converged
%       .iterations
%       .z
%       .A
%       .tofError

if nargin < 5
    longWay = false;
end

r1 = r1(:);
r2 = r2(:);

r1n = norm(r1);
r2n = norm(r2);

cosDtheta = dot(r1, r2) / (r1n * r2n);
cosDtheta = max(-1, min(1, cosDtheta));
dtheta = acos(cosDtheta);

cross12 = cross(r1, r2);
if longWay
    if cross12(3) >= 0
        dtheta = 2*pi - dtheta;
    end
else
    if cross12(3) < 0
        dtheta = 2*pi - dtheta;
    end
end

A = sin(dtheta) * sqrt(r1n * r2n / (1 - cos(dtheta)));

if abs(A) < 1e-14
    error('Lambert solver failed: geometry is degenerate.');
end

z = 0;
tol = 1e-10;
maxIter = 200;
converged = false;

for k = 1:maxIter
    [S, C] = stumpff(z);
    
    y = r1n + r2n + A * (z*S - 1) / sqrt(C);
    
    if A > 0 && y < 0
        z = z + 0.1;
        continue
    end
    
    x = sqrt(y / C);
    tof_z = (x^3 * S + A * sqrt(y)) / sqrt(mu);
    
    F = tof_z - tof;
    
    if abs(F) < tol
        converged = true;
        break
    end
    
    % Numerical derivative for robustness in first implementation
    dz = 1e-6;
    [Sp, Cp] = stumpff(z + dz);
    yp = r1n + r2n + A * ((z + dz)*Sp - 1) / sqrt(Cp);
    
    if A > 0 && yp < 0
        z = z + 0.1;
        continue
    end
    
    xp = sqrt(yp / Cp);
    tof_p = (xp^3 * Sp + A * sqrt(yp)) / sqrt(mu);
    
    dFdz = (tof_p - tof_z) / dz;
    
    zNew = z - F / dFdz;
    
    if ~isfinite(zNew)
        error('Lambert solver failed: iteration produced non-finite z.');
    end
    
    z = zNew;
end

if ~converged
    warning('Lambert solver did not fully converge.')
end

[S, C] = stumpff(z);
y = r1n + r2n + A * (z*S - 1) / sqrt(C);

f    = 1 - y / r1n;
g    = A * sqrt(y / mu);
gdot = 1 - y / r2n;

v1 = (r2 - f*r1) / g;
v2 = (gdot*r2 - r1) / g;

sol.v1 = v1;
sol.v2 = v2;
sol.converged = converged;
sol.iterations = k;
sol.z = z;
sol.A = A;
sol.tofError = F;

end

function [S, C] = stumpff(z)
%STUMPFF Computes Stumpff functions S(z) and C(z).

if z > 1e-12
    sz = sqrt(z);
    S = (sz - sin(sz)) / (sz^3);
    C = (1 - cos(sz)) / z;
elseif z < -1e-12
    sz = sqrt(-z);
    S = (sinh(sz) - sz) / (sz^3);
    C = (cosh(sz) - 1) / (-z);
else
    S = 1/6;
    C = 1/2;
end

end