function out = differentialCorrectionPlanarLyapunov(x0, vy0, mu, maxIter, tol)
%DIFFERENTIALCORRECTIONPLANARLYAPUNOV Correct planar Lyapunov seed.
%
% Initial condition assumed:
%   [x0; 0; 0; 0; vy0; 0]
%
% Symmetry condition at half-period crossing:
%   y = 0 crossing with xd = 0
%
% Corrects vy0 while keeping x0 fixed.
%
% Output fields:
%   .x0
%   .vy0
%   .state0
%   .halfPeriod
%   .period
%   .converged
%   .iterations
%   .stateHalf
%   .trajectory

if nargin < 4 || isempty(maxIter)
    maxIter = 20;
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
        error('No symmetry-plane crossing detected.');
    end

    % Use the first detected nontrivial crossing
    stateHalf = YE(1,1:6).';
    PhiHalf = reshape(YE(1,7:end), 6, 6);

    xdHalf = stateHalf(4);

    traj.t = tHist;
    traj.Y = YHist;

    if abs(xdHalf) < tol
        converged = true;
        break
    end

    % Sensitivity d(xd_half)/d(vy0)
    dxd_dvy0 = PhiHalf(4,5);

    if abs(dxd_dvy0) < 1e-12
        error('Differential correction became singular.');
    end

    vy = vy - xdHalf / dxd_dvy0;
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
% Detect nontrivial y=0 crossing after launch from symmetry plane
if t < 1e-6
    value = 1;   % avoid immediate trigger at t=0
else
    value = Y(2);
end
isterminal = 1;
direction = 1;
end