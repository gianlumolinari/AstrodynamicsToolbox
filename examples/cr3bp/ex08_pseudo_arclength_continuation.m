clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/AstroToolbox')
startup

% ==========================================================
% Pseudo-arclength continuation demo
% Earth-Moon L1 planar Lyapunov family
% ==========================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('Pseudo-Arclength Continuation Demo\n');
fprintf('============================================================\n');

% ----------------------------------------------------------
% Retrieve two nearby seeds from JPL
% ----------------------------------------------------------
data = astro.cr3bp.queryJPLPeriodicOrbits( ...
    'sys', 'earth-moon', ...
    'family', 'lyapunov', ...
    'libr', 1);

seed1 = astro.cr3bp.parseJPLPeriodicOrbit(data, 1);
seed2 = astro.cr3bp.parseJPLPeriodicOrbit(data, 2);

mu = seed1.mu;

fprintf('\nSeed 1:\n');
fprintf('  x0     = %.12f\n', seed1.state(1));
fprintf('  vy0    = %.12f\n', seed1.state(5));
fprintf('  Period = %.12f TU\n', seed1.period);

fprintf('\nSeed 2:\n');
fprintf('  x0     = %.12f\n', seed2.state(1));
fprintf('  vy0    = %.12f\n', seed2.state(5));
fprintf('  Period = %.12f TU\n', seed2.period);

% ----------------------------------------------------------
% Continuation settings
% ----------------------------------------------------------
nMembers = 500;
ds = 5e-4;

family = astro.cr3bp.continueFamilyPseudoArc(seed1, seed2, nMembers, ds, mu);

nPlot = numel(family);
if nPlot < 2
    error('Pseudo-arclength continuation produced too few family members.');
end

Cvals = [family.C];

fprintf('\nFamily summary:\n');
fprintf('  Stored members : %d\n', nPlot);
fprintf('  C range        : [%.12f, %.12f]\n', min(Cvals), max(Cvals));

L = astro.cr3bp.lagrangePoints(mu);

% ----------------------------------------------------------
% x-y family plot
% ----------------------------------------------------------
figure('Color','w');
hold on
axis equal
grid on
box on

cmap = parula(max(nPlot,16));

for k = 1:nPlot
    traj = family(k).traj;
    plot(traj(:,1), traj(:,2), 'Color', cmap(k,:), 'LineWidth', 1.2);
end

plot(-mu, 0, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
plot(1-mu, 0, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'k');
plot(L.L1(1), L.L1(2), 'rs', 'MarkerSize', 7, 'LineWidth', 1.4);

text(L.L1(1), L.L1(2), '  L1', 'FontSize', 10, 'FontWeight', 'bold');

colormap(cmap);
cb = colorbar;
ylabel(cb, 'Family index', 'FontSize', 12, 'FontWeight', 'bold');
clim([1 max(nPlot,2)]);

xlabel('x [-]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('y [-]', 'FontSize', 13, 'FontWeight', 'bold');
title('L1 Lyapunov Family (Pseudo-Arclength Continuation)', ...
    'FontSize', 15, 'FontWeight', 'bold');

% ----------------------------------------------------------
% Jacobi trend
% ----------------------------------------------------------
figure('Color','w');
hold on
grid on
box on

plot(1:nPlot, Cvals, 'o-', 'LineWidth', 1.2, 'MarkerSize', 6);

xlabel('Family member index', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Jacobi constant C', 'FontSize', 13, 'FontWeight', 'bold');
title('Jacobi Constant Along Pseudo-Arclength Family', ...
    'FontSize', 15, 'FontWeight', 'bold');

fprintf('\nDone.\n');