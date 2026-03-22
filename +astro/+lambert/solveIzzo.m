function sol = solveIzzo(r1, r2, tof, mu, longWay)
%SOLVEIZZO Single-revolution Izzo-style Lambert solver using Householder.
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

if (r1(1)*r2(2) - r1(2)*r2(1)) < 0
    lambda = -lambda;
    it1 = cross(ir1, ih);
    it2 = cross(ir2, ih);
else
    it1 = cross(ih, ir1);
    it2 = cross(ih, ir2);
end

if longWay
    lambda = -lambda;
    it1 = -it1;
    it2 = -it2;
end

% Non-dimensional target TOF
Tstar = sqrt(2*mu / s^3) * tof;

% Initial guess based on Izzo single-rev formulas
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
maxIter = 20;
converged = false;

for k = 1:maxIter
    tofData = astro.lambert.computeTOF_Izzo(x, lambda, 0);
    
    f   = tofData.T - Tstar;
    fp  = tofData.dTdx;
    fpp = tofData.d2Tdx2;
    fppp = tofData.d3Tdx3;
    
    if abs(f) < tol
        converged = true;
        break
    end
    
    % If too close to parabolic singular higher derivatives, fall back to Newton
    if ~isfinite(fpp) || ~isfinite(fppp)
        dx = -f / fp;
    else
        % Householder update as in Izzo paper discussion
        num = f * (fp^2 - 0.5*f*fpp);
        den = fp * (fp^2 - f*fpp) + (f^2 * fppp)/6;
        dx = -num / den;
    end
    
    xNew = x + dx;
    
    if ~isfinite(xNew)
        error('Izzo solver: non-finite iterate.');
    end
    
    x = xNew;
end

tofData = astro.lambert.computeTOF_Izzo(x, lambda, 0);
F = tofData.T - Tstar;
y = tofData.y;

% Velocity reconstruction
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
sol.backend = 'householder-analytical-derivatives';
end