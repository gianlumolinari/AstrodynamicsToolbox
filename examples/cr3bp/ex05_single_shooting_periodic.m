clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/AstroToolbox')
startup

% ==========================================================
% Single-shooting periodic correction demo
%
% Uses a JPL Earth-Moon L1 planar Lyapunov seed, perturbs it,
% and corrects it with a generic single-shooting periodicity method.
% ==========================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('Single-Shooting Periodic Correction Demo\n');
fprintf('============================================================\n');

% ----------------------------------------------------------
% Retrieve JPL seed
% ----------------------------------------------------------
data = astro.cr3bp.queryJPLPeriodicOrbits( ...
    'sys', 'earth-moon', ...
    'family', 'lyapunov', ...
    'libr', 1);

orbit = astro.cr3bp.parseJPLPeriodicOrbit(data, 1);

mu = orbit.mu;
x0True = orbit.state;
Ttrue = orbit.period;

fprintf('\nJPL seed:\n');
fprintf('  Period    : %.12f TU\n', Ttrue);
fprintf('  Jacobi    : %.12f\n', orbit.jacobi);
fprintf('  Stability : %.12f\n', orbit.stability);

% ----------------------------------------------------------
% Perturb the seed
% Keep x fixed as phase anchor
% ----------------------------------------------------------
x0Guess = x0True;
x0Guess(2) = x0Guess(2) + 1e-4;
x0Guess(4) = x0Guess(4) - 1e-4;
x0Guess(5) = x0Guess(5) + 1e-4;

Tguess = Ttrue * 1.0005;

freeIdx = [2 3 4 5 6];

out = astro.cr3bp.singleShootingCorrector( ...
    x0Guess, Tguess, mu, freeIdx, true, 20, 1e-11);

fprintf('\nSingle-shooting result:\n');
fprintf('  Converged  : %d\n', out.converged);
fprintf('  Iterations : %d\n', out.iterations);
fprintf('  Residual   : %.6e\n', out.residual);
fprintf('  Period     : %.12f TU\n', out.tf);

Ccorr = astro.cr3bp.jacobiConstant(out.x0.', mu);
fprintf('  Corrected C: %.12f\n', Ccorr);

% ----------------------------------------------------------
% Propagate corrected orbit
% ----------------------------------------------------------
opts.RelTol = 1e-12;
opts.AbsTol = 1e-12;
opts.Solver = 'ode113';

fullOut = astro.propagators.propagate( ...
    @(t,x) astro.cr3bp.cr3bpRHS(t, x, mu), ...
    [0, out.tf], out.x0, opts);

L = astro.cr3bp.lagrangePoints(mu);

figure('Color','w');
hold on
axis equal
grid on
box on

plot(fullOut.x(:,1), fullOut.x(:,2), 'b-', 'LineWidth', 1.5);
plot(-mu, 0, 'bo', 'MarkerSize', 10, 'LineWidth', 1.4, 'MarkerFaceColor', 'b');
plot(1-mu, 0, 'ko', 'MarkerSize', 7, 'LineWidth', 1.4, 'MarkerFaceColor', 'k');
plot(L.L1(1), L.L1(2), 'rs', 'MarkerSize', 7, 'LineWidth', 1.4);
plot(L.L2(1), L.L2(2), 'rs', 'MarkerSize', 7, 'LineWidth', 1.4);

text(L.L1(1), L.L1(2), '  L1', 'FontSize', 10, 'FontWeight', 'bold');
text(L.L2(1), L.L2(2), '  L2', 'FontSize', 10, 'FontWeight', 'bold');

xlabel('x [-]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('y [-]', 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('Single-Shooting Corrected Periodic Orbit (C = %.6f)', Ccorr), ...
    'FontSize', 15, 'FontWeight', 'bold');
legend('Orbit', 'Primary', 'Secondary', 'L1/L2', 'Location', 'best');

fprintf('\nDone.\n');