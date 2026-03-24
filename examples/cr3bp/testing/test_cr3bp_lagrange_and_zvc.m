clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/AstroToolbox')
startup

% ==========================================================
% CR3BP example:
%   - Earth-Moon system
%   - Lagrange points
%   - zero-velocity curves in the x-y plane
% ==========================================================

mu = 0.012150585609624;   % Earth-Moon CR3BP mass parameter

L = astro.cr3bp.lagrangePoints(mu);

fprintf('\n');
fprintf('============================================================\n');
fprintf('CR3BP Lagrange Points and Zero-Velocity Curves\n');
fprintf('============================================================\n');
fprintf('mu = %.15f\n', mu);

fprintf('\nLagrange points:\n');
fprintf('  L1 = [%.12f, %.12f, %.12f]\n', L.L1(1), L.L1(2), L.L1(3));
fprintf('  L2 = [%.12f, %.12f, %.12f]\n', L.L2(1), L.L2(2), L.L2(3));
fprintf('  L3 = [%.12f, %.12f, %.12f]\n', L.L3(1), L.L3(2), L.L3(3));
fprintf('  L4 = [%.12f, %.12f, %.12f]\n', L.L4(1), L.L4(2), L.L4(3));
fprintf('  L5 = [%.12f, %.12f, %.12f]\n', L.L5(1), L.L5(2), L.L5(3));

% Jacobi constants at the Lagrange points
CL1 = astro.cr3bp.jacobiConstant([L.L1.' 0 0 0], mu);
CL2 = astro.cr3bp.jacobiConstant([L.L2.' 0 0 0], mu);
CL3 = astro.cr3bp.jacobiConstant([L.L3.' 0 0 0], mu);
CL4 = astro.cr3bp.jacobiConstant([L.L4.' 0 0 0], mu);
CL5 = astro.cr3bp.jacobiConstant([L.L5.' 0 0 0], mu);

fprintf('\nJacobi constants at the Lagrange points:\n');
fprintf('  C(L1) = %.12f\n', CL1);
fprintf('  C(L2) = %.12f\n', CL2);
fprintf('  C(L3) = %.12f\n', CL3);
fprintf('  C(L4) = %.12f\n', CL4);
fprintf('  C(L5) = %.12f\n', CL5);

% ----------------------------------------------------------
% Zero-velocity curves
% ----------------------------------------------------------
xv = linspace(-1.5, 1.5, 800);
yv = linspace(-1.5, 1.5, 800);
[X, Y] = meshgrid(xv, yv);

Omega = astro.cr3bp.effectivePotential(X, Y, 0, mu);

% Choose a Jacobi level slightly below C(L1) to show neck opening
Cplot = CL1 - 0.02;
Z = 2*Omega - Cplot;

figure('Color','w');
hold on
axis equal
grid on
box on

% Forbidden/allowed boundary: Z = 0
contour(X, Y, Z, [0 0], 'LineWidth', 1.6);

% Primaries
plot(-mu, 0, 'bo', 'MarkerSize', 10, 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
plot(1-mu, 0, 'ko', 'MarkerSize', 7, 'LineWidth', 1.5, 'MarkerFaceColor', 'k');

% Lagrange points
plot(L.L1(1), L.L1(2), 'rs', 'MarkerSize', 7, 'LineWidth', 1.4);
plot(L.L2(1), L.L2(2), 'rs', 'MarkerSize', 7, 'LineWidth', 1.4);
plot(L.L3(1), L.L3(2), 'rs', 'MarkerSize', 7, 'LineWidth', 1.4);
plot(L.L4(1), L.L4(2), 'ms', 'MarkerSize', 7, 'LineWidth', 1.4);
plot(L.L5(1), L.L5(2), 'ms', 'MarkerSize', 7, 'LineWidth', 1.4);

text(L.L1(1), L.L1(2), '  L1', 'FontSize', 10, 'FontWeight', 'bold');
text(L.L2(1), L.L2(2), '  L2', 'FontSize', 10, 'FontWeight', 'bold');
text(L.L3(1), L.L3(2), '  L3', 'FontSize', 10, 'FontWeight', 'bold');
text(L.L4(1), L.L4(2), '  L4', 'FontSize', 10, 'FontWeight', 'bold');
text(L.L5(1), L.L5(2), '  L5', 'FontSize', 10, 'FontWeight', 'bold');

xlabel('x [-]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('y [-]', 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('Earth-Moon CR3BP Zero-Velocity Curve (C = %.4f)', Cplot), ...
    'FontSize', 15, 'FontWeight', 'bold');

legend('Zero-velocity curve', 'Earth', 'Moon', ...
       'L1/L2/L3', 'L4/L5', 'Location', 'best');

fprintf('\nDone.\n');