clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/AstroToolbox')
startup

% ==========================================================
% Multiple-shooting periodic correction demo
%
% Uses a JPL Earth-Moon L1 halo seed, samples it into nodes,
% perturbs the node states, and restores periodic continuity.
% ==========================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('Multiple-Shooting Periodic Correction Demo\n');
fprintf('============================================================\n');

% ----------------------------------------------------------
% Retrieve JPL halo seed
% ----------------------------------------------------------
data = astro.cr3bp.queryJPLPeriodicOrbits( ...
    'sys', 'earth-moon', ...
    'family', 'halo', ...
    'libr', 1, ...
    'branch', 'N');

orbit = astro.cr3bp.parseJPLPeriodicOrbit(data, 1);

mu = orbit.mu;
x0 = orbit.state;
T = orbit.period;

fprintf('\nJPL seed:\n');
fprintf('  Period    : %.12f TU\n', T);
fprintf('  Jacobi    : %.12f\n', orbit.jacobi);
fprintf('  Stability : %.12f\n', orbit.stability);

% ----------------------------------------------------------
% Sample seed orbit into nodes
% ----------------------------------------------------------
N = 8;
dtSeg = (T / N) * ones(1, N);

opts.RelTol = 1e-12;
opts.AbsTol = 1e-12;
opts.Solver = 'ode113';

seedOut = astro.propagators.propagate( ...
    @(t,x) astro.cr3bp.eomCR3BP(t, x, mu), ...
    linspace(0, T, N+1), x0, opts);

Xnodes = seedOut.x(1:N,:).';

% Keep node 1 fixed, perturb nodes 2..N
rng(1);
Xnodes(:,2:end) = Xnodes(:,2:end) + 1e-5 * randn(6, N-1);

% ----------------------------------------------------------
% Run multiple-shooting corrector
% ----------------------------------------------------------
out = astro.cr3bp.multipleShootingCorrector(Xnodes, dtSeg, mu, 20, 1e-10);

fprintf('\nMultiple-shooting result:\n');
fprintf('  Converged  : %d\n', out.converged);
fprintf('  Iterations : %d\n', out.iterations);
fprintf('  Residual   : %.6e\n', out.residual);

% ----------------------------------------------------------
% Plot corrected segments
% ----------------------------------------------------------
figure('Color','w');
hold on
grid on
box on
axis equal

for i = 1:numel(out.segments)
    X = out.segments{i}.x;
    plot3(X(:,1), X(:,2), X(:,3), 'b-', 'LineWidth', 1.3);
end

plot3(-mu, 0, 0, 'bo', 'MarkerSize', 10, 'LineWidth', 1.4, 'MarkerFaceColor', 'b');
plot3(1-mu, 0, 0, 'ko', 'MarkerSize', 7, 'LineWidth', 1.4, 'MarkerFaceColor', 'k');

xlabel('x [-]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('y [-]', 'FontSize', 13, 'FontWeight', 'bold');
zlabel('z [-]', 'FontSize', 13, 'FontWeight', 'bold');
title('Multiple-Shooting Corrected Halo Orbit', 'FontSize', 15, 'FontWeight', 'bold');
view(35, 25)

fprintf('\nDone.\n');