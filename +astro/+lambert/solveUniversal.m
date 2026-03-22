function sol = solveUniversal(r1, r2, tof, mu, longWay)
%SOLVEUNIVERSAL Robust 0-rev Lambert solver using universal variables.
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

% -------------------------------------------------------------------------
% Solve F(z) = TOF(z) - tof using bracketing + fzero
% -------------------------------------------------------------------------
Ffun = @(z) lambertTOFResidual(z, r1n, r2n, A, tof, mu);

% Initial bracket
zL = -4;
zU =  4;

FL = safeEval(Ffun, zL);
FU = safeEval(Ffun, zU);

maxExpand = 50;
expandCount = 0;

while ~(isfinite(FL) && isfinite(FU) && sign(FL) ~= sign(FU))
    expandCount = expandCount + 1;
    if expandCount > maxExpand
        error('Universal Lambert solver failed to bracket the root.');
    end
    
    zL = zL - 2;
    zU = zU + 2;
    
    FL = safeEval(Ffun, zL);
    FU = safeEval(Ffun, zU);
end

opts = optimset('TolX',1e-12,'Display','off');
[z, ~, exitflag, output] = fzero(Ffun, [zL zU], opts);

converged = (exitflag > 0);

% Final reconstruction
[S, C] = stumpff(z);
y = r1n + r2n + A * (z*S - 1) / sqrt(C);

f    = 1 - y / r1n;
g    = A * sqrt(y / mu);
gdot = 1 - y / r2n;

v1 = (r2 - f*r1) / g;
v2 = (gdot*r2 - r1) / g;

tofError = lambertTOFResidual(z, r1n, r2n, A, tof, mu);

sol.v1 = v1;
sol.v2 = v2;
sol.converged = converged;
sol.iterations = output.iterations;
sol.z = z;
sol.A = A;
sol.tofError = tofError;

end

% =========================================================================
function F = lambertTOFResidual(z, r1n, r2n, A, tofTarget, mu)
[S, C] = stumpff(z);

if C <= 0 || ~isfinite(C)
    F = NaN;
    return
end

y = r1n + r2n + A * (z*S - 1) / sqrt(C);

if y <= 0 || ~isfinite(y)
    F = NaN;
    return
end

x = sqrt(y / C);
tof = (x^3 * S + A * sqrt(y)) / sqrt(mu);

F = tof - tofTarget;
end

% =========================================================================
function val = safeEval(fun, z)
try
    val = fun(z);
catch
    val = NaN;
end
end

% =========================================================================
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