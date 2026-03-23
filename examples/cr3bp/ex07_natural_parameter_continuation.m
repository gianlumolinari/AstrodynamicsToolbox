clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/AstroToolbox')
startup

% ==========================================================
% Natural parameter continuation demo
% Earth-Moon L1 Lyapunov family
% ==========================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('Natural Parameter Continuation Demo\n');
fprintf('============================================================\n');

% ----------------------------------------------------------
% Get seed from JPL
% ----------------------------------------------------------
data = astro.cr3bp.queryJPLPeriodicOrbits( ...
    'sys','earth-moon', ...
    'family','lyapunov', ...
    'libr',1);

seed = astro.cr3bp.parseJPLPeriodicOrbit(data, 1);
mu = seed.mu;

fprintf('\nSeed orbit:\n');
fprintf('  Period = %.12f TU\n', seed.period);

% ----------------------------------------------------------
% Continuation settings
% ----------------------------------------------------------
nMembers = 30;
stepSize = -1e-3;
nNodes   = 10;

% ----------------------------------------------------------
% Run continuation
% ----------------------------------------------------------
family = astro.cr3bp.continueFamilyNatural( ...
    seed, nMembers, stepSize, nNodes, mu);

if isempty(family)
    error('No converged family members were generated.');
end

nPlot = numel(family);
Cvals = NaN(1,nPlot);

for k = 1:nPlot
    Cvals(k) = astro.cr3bp.jacobiConstant(family(k).x0.', mu);
end

fprintf('\nConverged family summary:\n');
fprintf('  Number of converged members : %d\n', nPlot);
fprintf('  Jacobi range               : [%.12f, %.12f]\n', min(Cvals), max(Cvals));

% ----------------------------------------------------------
% Plot family in x-y plane
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

L = astro.cr3bp.lagrangePoints(mu);

plot(-mu, 0, 'bo', 'MarkerSize', 10, 'MarkerFaceColor','b');
plot(1-mu, 0, 'ko', 'MarkerSize', 7, 'MarkerFaceColor','k');
plot(L.L1(1), L.L1(2), 'rs', 'MarkerSize', 7, 'LineWidth', 1.4);

text(L.L1(1), L.L1(2), '  L1', 'FontSize', 10, 'FontWeight', 'bold');

colormap(cmap);
cb = colorbar;
ylabel(cb, 'Family index', 'FontSize', 12, 'FontWeight', 'bold');
clim([1 max(nPlot,2)]);

xlabel('x [-]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('y [-]', 'FontSize', 13, 'FontWeight', 'bold');
title('L1 Lyapunov Family (Natural Continuation)', ...
    'FontSize', 15, 'FontWeight', 'bold');

% ----------------------------------------------------------
% Jacobi vs family member
% ----------------------------------------------------------
figure('Color','w');
hold on
grid on
box on

plot(1:nPlot, Cvals, 'o-', 'LineWidth', 1.2, 'MarkerSize', 6);

xlabel('Family member index', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Jacobi constant C', 'FontSize', 13, 'FontWeight', 'bold');
title('Jacobi Constant Along Natural Continuation Family', ...
    'FontSize', 15, 'FontWeight', 'bold');

fprintf('\nDone.\n');