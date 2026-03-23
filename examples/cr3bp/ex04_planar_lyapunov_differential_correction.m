clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/AstroToolbox')
startup

% ==========================================================
% Planar Lyapunov differential correction demo
% Earth-Moon CR3BP
% ==========================================================

mu = 0.012150585609624;
L = astro.cr3bp.lagrangePoints(mu);

% Better seed near L1
x0 = L.L1(1) - 0.005;
vy0 = 0.045;

out = astro.cr3bp.differentialCorrectionPlanarLyapunov(x0, vy0, mu, 30, 1e-10);

fprintf('\n');
fprintf('============================================================\n');
fprintf('Planar Lyapunov Differential Correction Demo\n');
fprintf('============================================================\n');
fprintf('Converged      : %d\n', out.converged);
fprintf('Iterations     : %d\n', out.iterations);
fprintf('Corrected vy0  : %.12f\n', out.vy0);
fprintf('Half period    : %.12f TU\n', out.halfPeriod);
fprintf('Full period    : %.12f TU\n', out.period);
fprintf('xd at crossing : %.6e\n', out.stateHalf(4));

opts.RelTol = 1e-12;
opts.AbsTol = 1e-12;
opts.Solver = 'ode113';

fullOut = astro.propagators.propagate( ...
    @(t,x) astro.cr3bp.eomCR3BP(t, x, mu), ...
    [0, out.period], out.state0, opts);

C = astro.cr3bp.jacobiConstant(out.state0.', mu);

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
title(sprintf('Planar Lyapunov Orbit near L1 (C = %.6f)', C), ...
    'FontSize', 15, 'FontWeight', 'bold');
legend('Orbit', 'Primary', 'Secondary', 'L1/L2', 'Location', 'best');

fprintf('\nDone.\n');