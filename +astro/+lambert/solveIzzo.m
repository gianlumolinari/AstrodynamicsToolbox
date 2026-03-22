function sol = solveIzzo(r1, r2, tof, mu, longWay)
%SOLVEIZZO Single-revolution Izzo-style Lambert solver.
%
% INPUTS
%   r1      : initial position vector [km]
%   r2      : final position vector [km]
%   tof     : time of flight [s]
%   mu      : gravitational parameter [km^3/s^2]
%   longWay : logical flag
%
% OUTPUT
%   sol : struct with fields
%       .v1
%       .v2
%       .converged
%       .iterations
%       .tofError
%       .x
%       .y
%       .lambda
%       .method
%       .backend

if nargin < 5
    longWay = false;
end

r1 = r1(:);
r2 = r2(:);

r1n = norm(r1);
r2n = norm(r2);

cvec = r2 - r1;
c = norm(cvec);
s = 0.5 * (r1n + r2n + c);

ir1 = r1 / r1n;
ir2 = r2 / r2n;

ih = cross(ir1, ir2);
ihNorm = norm(ih);

if ihNorm < 1e-14
    error('Izzo solver: degenerate geometry (collinear vectors).');
end

ih = ih / ihNorm;

lambda2 = 1 - c/s;
lambda = sqrt(lambda2);

% Orientation / transfer direction
if (r1(1)*r2(2) - r1(2)*r2(1)) < 0
    lambda = -lambda;
    it1 = cross(ir1, ih);
    it2 = cross(ir2, ih);
else
    it1 = cross(ih, ir1);
    it2 = cross(ih, ir2);
end

% Apply long-way choice by flipping lambda sign
if longWay
    lambda = -lambda;
    it1 = -it1;
    it2 = -it2;
end

% Non-dimensional TOF
Tstar = sqrt(2*mu / s^3) * tof;

% Single-revolution initial guess from Izzo Eq. (30)-type logic
T0 = acos(lambda) + lambda*sqrt(1 - lambda^2);
T1 = (2/3) * (1 - lambda^3);

if Tstar >= T0
    x = (T0 / Tstar)^(2/3) - 1;
elseif Tstar < T1
    x = 2.5 * T1 * (T1 - Tstar) / (Tstar * (1 - lambda^5)) + 1;
else
    x = (T0 / Tstar)^(log(2) / log(T1 / T0)) - 1;
end

tol = 1e-12;
maxIter = 50;
converged = false;

for k = 1:maxIter
    [T, y] = astro.lambert.computeTOF_Izzo(x, lambda, 0);
    F = T - Tstar;

    if abs(F) < tol
        converged = true;
        break
    end

    % Numerical derivative for now
    dx = 1e-7;
    [Tp, ~] = astro.lambert.computeTOF_Izzo(x + dx, lambda, 0);
    dFdx = (Tp - T) / dx;

    if abs(dFdx) < 1e-14 || ~isfinite(dFdx)
        error('Izzo solver: derivative breakdown.');
    end

    xNew = x - F / dFdx;

    if ~isfinite(xNew)
        error('Izzo solver: non-finite iterate.');
    end

    x = xNew;
end

[T, y] = astro.lambert.computeTOF_Izzo(x, lambda, 0);
F = T - Tstar;

% Velocity reconstruction from Izzo algorithm
gamma = sqrt(mu * s / 2);
rho   = (r1n - r2n) / c;
sigma = sqrt(1 - rho^2);

Vr1 = gamma * ((lambda*y - x) - rho*(lambda*y + x)) / r1n;
Vr2 = -gamma * ((lambda*y - x) + rho*(lambda*y + x)) / r2n;
Vt1 = gamma * sigma * (y + lambda*x) / r1n;
Vt2 = gamma * sigma * (y + lambda*x) / r2n;

v1 = Vr1*ir1 + Vt1*it1;
v2 = Vr2*ir2 + Vt2*it2;

sol.v1 = v1;
sol.v2 = v2;
sol.converged = converged;
sol.iterations = k;
sol.tofError = F;
sol.x = x;
sol.y = y;
sol.lambda = lambda;
sol.method = 'Izzo';
sol.backend = 'real-x-iteration-corrected';
end