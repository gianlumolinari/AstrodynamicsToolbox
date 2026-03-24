clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/AstroToolbox')
startup

% ==========================================================
% Manifold tube visualization demo
% Earth-Moon L1 planar Lyapunov family
%
% Goal:
%   Produce a cleaner, denser visualization of stable and
%   unstable manifold bundles ("tubes") from one periodic orbit.
% ==========================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('Manifold Tube Visualization Demo\n');
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
% Dense manifold seed set
% ----------------------------------------------------------
nPhase = 80;          % denser than before
epsManifold = 5e-7;   % slightly smaller perturbation

seeds = astro.cr3bp.manifoldSeeds(x0, T, mu, nPhase, epsManifold);

fprintf('\nGenerated %d phase points for manifold seeding.\n', nPhase);

% ----------------------------------------------------------
% Propagation settings
% ----------------------------------------------------------
tUnstable = 8.0;
tStable   = -8.0;

unstablePlus  = cell(nPhase,1);
unstableMinus = cell(nPhase,1);
stablePlus    = cell(nPhase,1);
stableMinus   = cell(nPhase,1);

fprintf('\nPropagating manifold bundles...\n');
for k = 1:nPhase
    unstablePlus{k}  = astro.cr3bp.propagateManifold(seeds(k).xUnstablePlus,  tUnstable, mu);
    unstableMinus{k} = astro.cr3bp.propagateManifold(seeds(k).xUnstableMinus, tUnstable, mu);

    stablePlus{k}    = astro.cr3bp.propagateManifold(seeds(k).xStablePlus,    tStable,   mu);
    stableMinus{k}   = astro.cr3bp.propagateManifold(seeds(k).xStableMinus,   tStable,   mu);
end

L = astro.cr3bp.lagrangePoints(mu);

% ----------------------------------------------------------
% Plot 1: global manifold tubes in x-y plane
% ----------------------------------------------------------
figure('Color','w');
hold on
axis equal
grid on
box on

% Periodic orbit
plot(family(kMember).traj(:,1), family(kMember).traj(:,2), ...
    'k-', 'LineWidth', 2.0);

% Stable bundles (greens)
for k = 1:nPhase
    X = stablePlus{k}.x;
    plot(X(:,1), X(:,2), '-', 'Color', [0.10 0.60 0.10], 'LineWidth', 0.8);

    X = stableMinus{k}.x;
    plot(X(:,1), X(:,2), '-', 'Color', [0.10 0.82 0.10], 'LineWidth', 0.8);
end

% Unstable bundles (reds)
for k = 1:nPhase
    X = unstablePlus{k}.x;
    plot(X(:,1), X(:,2), '-', 'Color', [0.80 0.10 0.10], 'LineWidth', 0.8);

    X = unstableMinus{k}.x;
    plot(X(:,1), X(:,2), '-', 'Color', [0.95 0.45 0.45], 'LineWidth', 0.8);
end

% Primaries and L1
plot(-mu, 0, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
plot(1-mu, 0, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'k');
plot(L.L1(1), L.L1(2), 'rs', 'MarkerSize', 7, 'LineWidth', 1.4);

text(L.L1(1), L.L1(2), '  L1', 'FontSize', 10, 'FontWeight', 'bold');

xlabel('$x~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('$y~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('Stable and Unstable Manifold Tubes from an $L_1$ Lyapunov Orbit ($C = %.6f$)', C), ...
    'Interpreter', 'latex', 'FontSize', 15, 'FontWeight', 'bold');

legend('Periodic orbit', ...
       'Stable manifold (+)', 'Stable manifold (-)', ...
       'Unstable manifold (+)', 'Unstable manifold (-)', ...
       'Primary', 'Secondary', '$L_1$', ...
       'Interpreter', 'latex', 'Location', 'best');

% ----------------------------------------------------------
% Plot 2: zoom near L1 / orbit neighborhood
% ----------------------------------------------------------
figure('Color','w');
hold on
axis equal
grid on
box on

plot(family(kMember).traj(:,1), family(kMember).traj(:,2), ...
    'k-', 'LineWidth', 2.0);

for k = 1:nPhase
    X = stablePlus{k}.x;
    plot(X(:,1), X(:,2), '-', 'Color', [0.10 0.60 0.10], 'LineWidth', 0.8);

    X = stableMinus{k}.x;
    plot(X(:,1), X(:,2), '-', 'Color', [0.10 0.82 0.10], 'LineWidth', 0.8);

    X = unstablePlus{k}.x;
    plot(X(:,1), X(:,2), '-', 'Color', [0.80 0.10 0.10], 'LineWidth', 0.8);

    X = unstableMinus{k}.x;
    plot(X(:,1), X(:,2), '-', 'Color', [0.95 0.45 0.45], 'LineWidth', 0.8);
end

plot(-mu, 0, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
plot(1-mu, 0, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'k');
plot(L.L1(1), L.L1(2), 'rs', 'MarkerSize', 7, 'LineWidth', 1.4);

xlim([L.L1(1)-0.08, L.L1(1)+0.08]);
ylim([-0.08, 0.08]);

xlabel('$x~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('$y~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
title('Local View of the Manifold Tube Structure Near $L_1$', ...
    'Interpreter', 'latex', 'FontSize', 15, 'FontWeight', 'bold');

% ----------------------------------------------------------
% Plot 3: state-space style view (x, y, ydot)
% This often makes the "tube" appearance clearer even for
% planar Lyapunov manifolds.
% ----------------------------------------------------------
figure('Color','w');
hold on
grid on
box on

% periodic orbit in reduced state space
plot3(family(kMember).traj(:,1), family(kMember).traj(:,2), family(kMember).traj(:,5), ...
    'k-', 'LineWidth', 2.0);

for k = 1:nPhase
    X = stablePlus{k}.x;
    plot3(X(:,1), X(:,2), X(:,5), '-', 'Color', [0.10 0.60 0.10], 'LineWidth', 0.8);

    X = stableMinus{k}.x;
    plot3(X(:,1), X(:,2), X(:,5), '-', 'Color', [0.10 0.82 0.10], 'LineWidth', 0.8);

    X = unstablePlus{k}.x;
    plot3(X(:,1), X(:,2), X(:,5), '-', 'Color', [0.80 0.10 0.10], 'LineWidth', 0.8);

    X = unstableMinus{k}.x;
    plot3(X(:,1), X(:,2), X(:,5), '-', 'Color', [0.95 0.45 0.45], 'LineWidth', 0.8);
end

xlabel('$x~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('$y~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
zlabel('$\dot{y}~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
title('Manifold Bundles in Reduced State Space $(x,y,\dot{y})$', ...
    'Interpreter', 'latex', 'FontSize', 15, 'FontWeight', 'bold');
view(35, 25)

fprintf('\nDone.\n');
