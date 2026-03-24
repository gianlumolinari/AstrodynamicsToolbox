clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/AstroToolbox')
startup

% ==========================================================
% Poincare section demo for Lyapunov-orbit manifolds
% Earth-Moon L1 planar Lyapunov family
% ==========================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('Poincare Section Demo\n');
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
nPhase = 20;
epsManifold = 1e-6;

seeds = astro.cr3bp.manifoldSeeds(x0, T, mu, nPhase, epsManifold);

% Build seed matrices
XunstablePlus  = zeros(6, nPhase);
XunstableMinus = zeros(6, nPhase);
XstablePlus    = zeros(6, nPhase);
XstableMinus   = zeros(6, nPhase);

for k = 1:nPhase
    XunstablePlus(:,k)  = seeds(k).xUnstablePlus;
    XunstableMinus(:,k) = seeds(k).xUnstableMinus;
    XstablePlus(:,k)    = seeds(k).xStablePlus;
    XstableMinus(:,k)   = seeds(k).xStableMinus;
end

fprintf('\nGenerated %d phase points for manifold seeding.\n', nPhase);

% ----------------------------------------------------------
% Define section plane
% ----------------------------------------------------------
xSection = 1.0 - mu;   % section near Moon's x-location
tUnstable = 6.0;
tStable   = -6.0;

dirUnstable = 0;
dirStable   = 0;

% ----------------------------------------------------------
% Collect section points
% ----------------------------------------------------------
secUplus  = astro.cr3bp.collectPoincareSection(XunstablePlus,  tUnstable, mu, xSection, dirUnstable);
secUminus = astro.cr3bp.collectPoincareSection(XunstableMinus, tUnstable, mu, xSection, dirUnstable);
secSplus  = astro.cr3bp.collectPoincareSection(XstablePlus,    tStable,   mu, xSection, dirStable);
secSminus = astro.cr3bp.collectPoincareSection(XstableMinus,   tStable,   mu, xSection, dirStable);

fprintf('\nSection hits:\n');
fprintf('  Unstable (+) : %d\n', size(secUplus.points,2));
fprintf('  Unstable (-) : %d\n', size(secUminus.points,2));
fprintf('  Stable   (+) : %d\n', size(secSplus.points,2));
fprintf('  Stable   (-) : %d\n', size(secSminus.points,2));

% ----------------------------------------------------------
% Plot configuration-space manifolds and section plane
% ----------------------------------------------------------
L = astro.cr3bp.lagrangePoints(mu);

figure('Color','w');
hold on
axis equal
grid on
box on

plot(family(kMember).traj(:,1), family(kMember).traj(:,2), 'k-', 'LineWidth', 1.8);

for k = 1:nPhase
    out = secUplus.trajectories{k};
    plot(out.x(:,1), out.x(:,2), '-', 'Color', [0.8 0.1 0.1], 'LineWidth', 1.0);

    out = secUminus.trajectories{k};
    plot(out.x(:,1), out.x(:,2), '-', 'Color', [0.95 0.4 0.4], 'LineWidth', 1.0);

    out = secSplus.trajectories{k};
    plot(out.x(:,1), out.x(:,2), '-', 'Color', [0.1 0.6 0.1], 'LineWidth', 1.0);

    out = secSminus.trajectories{k};
    plot(out.x(:,1), out.x(:,2), '-', 'Color', [0.1 0.8 0.1], 'LineWidth', 1.0);
end

xline(xSection, '--', 'LineWidth', 1.2);

plot(-mu, 0, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
plot(1-mu, 0, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'k');
plot(L.L1(1), L.L1(2), 'rs', 'MarkerSize', 7, 'LineWidth', 1.4);

xlabel('$x~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('$y~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('Manifolds and Section $x = %.6f$', xSection), ...
    'Interpreter', 'latex', 'FontSize', 15, 'FontWeight', 'bold');

legend('Periodic orbit', ...
       'Unstable (+)', 'Unstable (-)', ...
       'Stable (+)', 'Stable (-)', ...
       'Section plane', 'Primary', 'Secondary', 'L1', ...
       'Location', 'best');

% ----------------------------------------------------------
% Plot Poincare section in (y, yd)
% ----------------------------------------------------------
figure('Color','w');
hold on
grid on
box on

if ~isempty(secUplus.points)
    plot(secUplus.points(2,:), secUplus.points(5,:), 'ro', 'MarkerSize', 6, 'LineWidth', 1.2);
end
if ~isempty(secUminus.points)
    plot(secUminus.points(2,:), secUminus.points(5,:), 'rs', 'MarkerSize', 6, 'LineWidth', 1.2);
end
if ~isempty(secSplus.points)
    plot(secSplus.points(2,:), secSplus.points(5,:), 'go', 'MarkerSize', 6, 'LineWidth', 1.2);
end
if ~isempty(secSminus.points)
    plot(secSminus.points(2,:), secSminus.points(5,:), 'gs', 'MarkerSize', 6, 'LineWidth', 1.2);
end

xlabel('$y$ at section $[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('$\dot{y}$ at section $[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('Poincaré Section at $x = %.6f$', xSection), ...
    'Interpreter', 'latex', 'FontSize', 15, 'FontWeight', 'bold');
legend('Unstable (+)', 'Unstable (-)', 'Stable (+)', 'Stable (-)', ...
    'Location', 'best');

% ----------------------------------------------------------
% Optional second section plot in (y, xd)
% ----------------------------------------------------------
figure('Color','w');
hold on
grid on
box on

if ~isempty(secUplus.points)
    plot(secUplus.points(2,:), secUplus.points(4,:), 'ro', 'MarkerSize', 6, 'LineWidth', 1.2);
end
if ~isempty(secUminus.points)
    plot(secUminus.points(2,:), secUminus.points(4,:), 'rs', 'MarkerSize', 6, 'LineWidth', 1.2);
end
if ~isempty(secSplus.points)
    plot(secSplus.points(2,:), secSplus.points(4,:), 'go', 'MarkerSize', 6, 'LineWidth', 1.2);
end
if ~isempty(secSminus.points)
    plot(secSminus.points(2,:), secSminus.points(4,:), 'gs', 'MarkerSize', 6, 'LineWidth', 1.2);
end

xlabel('$y$ at section $[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('$\dot{x}$ at section $[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('Alternative Poincaré Section at $x = %.6f$', xSection), ...
    'Interpreter', 'latex', 'FontSize', 15, 'FontWeight', 'bold');
legend('Unstable (+)', 'Unstable (-)', 'Stable (+)', 'Stable (-)', ...
    'Location', 'best');

fprintf('\nDone.\n');