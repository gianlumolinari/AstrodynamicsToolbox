clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/AstroToolbox')
startup

% ==========================================================
% Halo family pseudo-arclength continuation demo
% Earth-Moon L1 northern halo family
% ==========================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('Halo Family Continuation Demo\n');
fprintf('============================================================\n');

% ----------------------------------------------------------
% Retrieve two nearby JPL L1 northern halo seeds
% ----------------------------------------------------------
data = astro.cr3bp.queryJPLPeriodicOrbits( ...
    'sys', 'earth-moon', ...
    'family', 'halo', ...
    'libr', 1, ...
    'branch', 'N');

seed1 = astro.cr3bp.parseJPLPeriodicOrbit(data, 1);
seed2 = astro.cr3bp.parseJPLPeriodicOrbit(data, 2);

mu = seed1.mu;

fprintf('\nSeed 1:\n');
fprintf('  x0     = %.12f\n', seed1.state(1));
fprintf('  z0     = %.12f\n', seed1.state(3));
fprintf('  vy0    = %.12f\n', seed1.state(5));
fprintf('  Period = %.12f TU\n', seed1.period);

fprintf('\nSeed 2:\n');
fprintf('  x0     = %.12f\n', seed2.state(1));
fprintf('  z0     = %.12f\n', seed2.state(3));
fprintf('  vy0    = %.12f\n', seed2.state(5));
fprintf('  Period = %.12f TU\n', seed2.period);

% ----------------------------------------------------------
% Continue family
% ----------------------------------------------------------
nMembers = 30;
ds = 5e-4;

family = astro.cr3bp.continueHaloPseudoArc(seed1, seed2, nMembers, ds, mu);

nFam = numel(family);
if nFam < 2
    error('Too few halo family members generated.');
end

Cvals = [family.C];
Tvals = [family.period];
z0vals = arrayfun(@(s) s.u(2), family);

fprintf('\nHalo family summary:\n');
fprintf('  Stored members : %d\n', nFam);
fprintf('  C range        : [%.12f, %.12f]\n', min(Cvals), max(Cvals));
fprintf('  T range [TU]   : [%.12f, %.12f]\n', min(Tvals), max(Tvals));
fprintf('  z0 range       : [%.12f, %.12f]\n', min(z0vals), max(z0vals));

L = astro.cr3bp.lagrangePoints(mu);

% ----------------------------------------------------------
% 3D family plot
% ----------------------------------------------------------
figure('Color','w');
hold on
grid on
box on
axis equal

cmap = parula(max(nFam,16));

for k = 1:nFam
    X = family(k).traj;
    plot3(X(:,1), X(:,2), X(:,3), 'Color', cmap(k,:), 'LineWidth', 1.2);
end

plot3(-mu, 0, 0, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
plot3(1-mu, 0, 0, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'k');
plot3(L.L1(1), L.L1(2), 0, 'rs', 'MarkerSize', 7, 'LineWidth', 1.4);

xlabel('$x~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('$y~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
zlabel('$z~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
title('$L_1$ Northern Halo Family (Pseudo-Arclength Continuation)', ...
    'Interpreter', 'latex', 'FontSize', 15, 'FontWeight', 'bold');
view(35, 25)

colormap(cmap);
cb = colorbar;
ylabel(cb, 'Family index', 'FontSize', 12, 'FontWeight', 'bold');
clim([1 max(nFam,2)]);

% ----------------------------------------------------------
% x-z projection
% ----------------------------------------------------------
figure('Color','w');
hold on
grid on
box on
axis equal

for k = 1:nFam
    X = family(k).traj;
    plot(X(:,1), X(:,3), 'Color', cmap(k,:), 'LineWidth', 1.2);
end

plot(-mu, 0, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
plot(1-mu, 0, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'k');
plot(L.L1(1), 0, 'rs', 'MarkerSize', 7, 'LineWidth', 1.4);

xlabel('$x~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('$z~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
title('$x$-$z$ Projection of the $L_1$ Halo Family', ...
    'Interpreter', 'latex', 'FontSize', 15, 'FontWeight', 'bold');

% ----------------------------------------------------------
% Jacobi and period trends
% ----------------------------------------------------------
figure('Color','w');

subplot(1,2,1)
hold on
grid on
box on
plot(1:nFam, Cvals, 'o-', 'LineWidth', 1.2, 'MarkerSize', 6);
xlabel('Family member index', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Jacobi constant $C$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
title('$C$ Along Halo Family', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');

subplot(1,2,2)
hold on
grid on
box on
plot(1:nFam, Tvals, 'o-', 'LineWidth', 1.2, 'MarkerSize', 6);
xlabel('Family member index', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Period [TU]', 'FontSize', 12, 'FontWeight', 'bold');
title('Period Along Halo Family', 'Interpreter', 'latex','FontSize', 13, 'FontWeight', 'bold');

fprintf('\nDone.\n');