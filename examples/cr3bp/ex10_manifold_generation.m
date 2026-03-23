clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/AstroToolbox')
startup

% ==========================================================
% Stable and unstable manifold generation demo
% Earth-Moon L1 planar Lyapunov family
% ==========================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('Manifold Generation Demo\n');
fprintf('============================================================\n');

% ----------------------------------------------------------
% Build a small pseudo-arclength family and pick one member
% ----------------------------------------------------------
data = astro.cr3bp.queryJPLPeriodicOrbits( ...
    'sys', 'earth-moon', ...
    'family', 'lyapunov', ...
    'libr', 1);

seed1 = astro.cr3bp.parseJPLPeriodicOrbit(data, 1);
seed2 = astro.cr3bp.parseJPLPeriodicOrbit(data, 2);

mu = seed1.mu;

family = astro.cr3bp.continueFamilyPseudoArc(seed1, seed2, 20, 5e-4, mu);

kMember = 10;
if kMember > numel(family)
    error('Requested family member exceeds generated family size.');
end

x0 = family(kMember).state0(:);
T  = family(kMember).period;
C  = family(kMember).C;

fprintf('\nSelected family member %d\n', kMember);
fprintf('  Period = %.12f TU\n', T);
fprintf('  Jacobi = %.12f\n', C);

% ----------------------------------------------------------
% Generate manifold seeds
% ----------------------------------------------------------
nPhase = 16;
epsManifold = 1e-6;

seeds = astro.cr3bp.manifoldSeeds(x0, T, mu, nPhase, epsManifold);

fprintf('\nGenerated %d phase points for manifold seeding.\n', numel(seeds));

% ----------------------------------------------------------
% Propagate manifolds
% ----------------------------------------------------------
tUnstable = 4.0;
tStable   = -4.0;

unstablePlus = cell(nPhase,1);
unstableMinus = cell(nPhase,1);
stablePlus = cell(nPhase,1);
stableMinus = cell(nPhase,1);

for k = 1:nPhase
    unstablePlus{k} = astro.cr3bp.propagateManifold(seeds(k).xUnstablePlus, tUnstable, mu);
    unstableMinus{k} = astro.cr3bp.propagateManifold(seeds(k).xUnstableMinus, tUnstable, mu);

    stablePlus{k} = astro.cr3bp.propagateManifold(seeds(k).xStablePlus, tStable, mu);
    stableMinus{k} = astro.cr3bp.propagateManifold(seeds(k).xStableMinus, tStable, mu);
end

% ----------------------------------------------------------
% Plot periodic orbit + manifolds
% ----------------------------------------------------------
L = astro.cr3bp.lagrangePoints(mu);

figure('Color','w');
hold on
axis equal
grid on
box on

% Plot the reference periodic orbit
plot(family(kMember).traj(:,1), family(kMember).traj(:,2), 'k-', 'LineWidth', 1.8);

% Stable manifolds
for k = 1:nPhase
    plot(stablePlus{k}.x(:,1), stablePlus{k}.x(:,2), '-', 'LineWidth', 1.0, ...
        'Color', [0.1 0.6 0.1]);
    plot(stableMinus{k}.x(:,1), stableMinus{k}.x(:,2), '-', 'LineWidth', 1.0, ...
        'Color', [0.1 0.8 0.1]);
end

% Unstable manifolds
for k = 1:nPhase
    plot(unstablePlus{k}.x(:,1), unstablePlus{k}.x(:,2), '-', 'LineWidth', 1.0, ...
        'Color', [0.8 0.1 0.1]);
    plot(unstableMinus{k}.x(:,1), unstableMinus{k}.x(:,2), '-', 'LineWidth', 1.0, ...
        'Color', [0.95 0.4 0.4]);
end

% Primaries and L1
plot(-mu, 0, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
plot(1-mu, 0, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'k');
plot(L.L1(1), L.L1(2), 'rs', 'MarkerSize', 7, 'LineWidth', 1.4);

text(L.L1(1), L.L1(2), '  L1', 'FontSize', 10, 'FontWeight', 'bold');

xlabel('x [-]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('y [-]', 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('Stable and Unstable Manifolds from L1 Lyapunov Orbit (C = %.6f)', C), ...
    'FontSize', 15, 'FontWeight', 'bold');

legend('Periodic orbit', ...
       'Stable manifold (+)', 'Stable manifold (-)', ...
       'Unstable manifold (+)', 'Unstable manifold (-)', ...
       'Primary', 'Secondary', 'L1', ...
       'Location', 'best');

fprintf('\nDone.\n');