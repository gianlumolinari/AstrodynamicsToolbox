function out = differentialCorrectionPlanarLyapunov(x0, vy0, mu, maxIter, tol)
%DIFFERENTIALCORRECTIONPLANARLYAPUNOV Differential correction for a planar
% Lyapunov orbit in the CR3BP.
%
% Initial state form:
%   [x0; 0; 0; 0; vy0; 0]
%
% Corrected variable:
%   vy0
%
% Symmetry condition at half period:
%   y = 0 crossing with xd = 0
%
% The correction includes the event-time dependence:
%   d(xd_f)/d(vy0) = Phi(4,5) - (xdd_f / yd_f)*Phi(2,5)

if nargin < 4 || isempty(maxIter)
    maxIter = 30;
end
if nargin < 5 || isempty(tol)
    tol = 1e-10;
end

vy = vy0;
converged = false;
traj = struct('t', [], 'Y', []);

for k = 1:maxIter
    Phi0 = eye(6);
    Y0 = [x0; 0; 0; 0; vy; 0; Phi0(:)];

    opts = odeset( ...
        'RelTol', 1e-12, ...
        'AbsTol', 1e-12, ...
        'Events', @(t,Y) localYCrossingEvent(t,Y));

    [tHist, YHist, tE, YE, ~] = ode113( ...
        @(t,Y) astro.cr3bp.variationalEOM(t, Y, mu), ...
        [0, 20], Y0, opts);

    if isempty(tE)
        error('No nontrivial y=0 crossing detected.');
    end

    stateHalf = YE(1,1:6).';
    PhiHalf = reshape(YE(1,7:end), 6, 6);

    xdf = stateHalf(4);
    ydf = stateHalf(5);

    traj.t = tHist;
    traj.Y = YHist;

    fprintf('    DC iter %2d: xdot(half) = %+ .3e, vy = %.16f, thalf = %.16f\n', ...
        k, xdf, vy, tE(1));

    if abs(xdf) < tol
        converged = true;
        break
    end

    % Accelerations at the crossing
    fHalf = astro.cr3bp.eomCR3BP(0, stateHalf, mu);
    xddf = fHalf(4);

    % Time-corrected sensitivity
    if abs(ydf) < 1e-12
        error('Differential correction became singular because ydot at crossing is too small.');
    end

    denom = PhiHalf(4,5) - (xddf / ydf) * PhiHalf(2,5);

    if abs(denom) < 1e-12
        error('Differential correction became singular because denominator is too small.');
    end

    dvy = -xdf / denom;
    vy = vy + dvy;
end

halfPeriod = tE(1);
period = 2 * halfPeriod;

out.x0 = x0;
out.vy0 = vy;
out.state0 = [x0; 0; 0; 0; vy; 0];
out.halfPeriod = halfPeriod;
out.period = period;
out.converged = converged;
out.iterations = k;
out.stateHalf = stateHalf;
out.trajectory = traj;

end

function [value, isterminal, direction] = localYCrossingEvent(t, Y)
% We start at y=0 with vy>0, so the half-period crossing is the NEXT
% y=0 crossing with negative direction (downward crossing).

if t < 1e-6
    value = 1;
else
    value = Y(2);
end

isterminal = 1;
direction = -1;
end