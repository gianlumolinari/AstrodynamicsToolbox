clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/AstroToolbox')
startup

% ==========================================================
% Halo manifold tube visualization demo
% Earth-Moon L1 northern halo family
% ==========================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('Halo Manifold Tube Visualization Demo\n');
fprintf('============================================================\n');

% ----------------------------------------------------------
% Build halo family and pick one member
% ----------------------------------------------------------
data = astro.cr3bp.queryJPLPeriodicOrbits( ...
    'sys', 'earth-moon', ...
    'family', 'halo', ...
    'libr', 1, ...
    'branch', 'N');

seed1 = astro.cr3bp.parseJPLPeriodicOrbit(data, 1);
seed2 = astro.cr3bp.parseJPLPeriodicOrbit(data, 2);

mu = seed1.mu;

family = astro.cr3bp.continueHaloPseudoArc(seed1, seed2, 20, 5e-4, mu);

kMember = 10;
if kMember > numel(family)
    error('Requested halo family member exceeds generated family size.');
end

x0 = family(kMember).state0(:);
T  = family(kMember).period;
C  = family(kMember).C;

fprintf('\nSelected halo family member %d\n', kMember);
fprintf('  Period = %.12f TU\n', T);
fprintf('  Jacobi = %.12f\n', C);

% ----------------------------------------------------------
% Tube seed settings
% ----------------------------------------------------------
nPhase = 36;
nCircle = 24;
epsManifold = 1e-7;

seeds = astro.manifolds.makeManifoldSeeds(x0, T, mu, nPhase, nCircle, epsManifold);

fprintf('\nGenerated halo tube seeds:\n');
fprintf('  Phase points : %d\n', nPhase);
fprintf('  Ring points  : %d\n', nCircle);

% ----------------------------------------------------------
% Propagation settings
% ----------------------------------------------------------
tUnstable = 2.5;
tStable   = -2.5;

unstableTraj = cell(nPhase, nCircle);
stableTraj   = cell(nPhase, nCircle);

fprintf('\nPropagating halo manifold tubes...\n');
for k = 1:nPhase
    for j = 1:nCircle
        unstableTraj{k,j} = astro.cr3bp.propagateManifold( ...
            seeds(k).xUnstableRing(:,j), tUnstable, mu);

        stableTraj{k,j} = astro.cr3bp.propagateManifold( ...
            seeds(k).xStableRing(:,j), tStable, mu);
    end
end

L = astro.cr3bp.lagrangePoints(mu);

% ----------------------------------------------------------
% Plot 1: 3D halo and manifold tubes
% ----------------------------------------------------------
figure('Color','w');
hold on
grid on
box on
axis equal

plot3(family(kMember).traj(:,1), family(kMember).traj(:,2), family(kMember).traj(:,3), ...
    'k-', 'LineWidth', 2.0);

for k = 1:nPhase
    for j = 1:nCircle
        X = unstableTraj{k,j}.x;
        plot3(X(:,1), X(:,2), X(:,3), '-', ...
            'Color', [0.88 0.25 0.25], 'LineWidth', 0.45);

        X = stableTraj{k,j}.x;
        plot3(X(:,1), X(:,2), X(:,3), '-', ...
            'Color', [0.25 0.72 0.40], 'LineWidth', 0.45);
    end
end

plot3(-mu, 0, 0, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
plot3(1-mu, 0, 0, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'k');
plot3(L.L1(1), 0, 0, 'rs', 'MarkerSize', 7, 'LineWidth', 1.4);

xlabel('$x~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('$y~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
zlabel('$z~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('Halo Orbit with Stable/Unstable Manifold Tubes $(C = %.6f)$', C), ...
    'Interpreter', 'latex', 'FontSize', 15, 'FontWeight', 'bold');
legend('Halo orbit', 'Unstable tube', 'Stable tube', 'Primary', 'Secondary', '$L_1$', ...
    'Interpreter', 'latex', 'Location', 'best');
view(35, 25)

% ----------------------------------------------------------
% Plot 2: x-z projection (often the nicest)
% ----------------------------------------------------------
figure('Color','w');
hold on
grid on
box on
axis equal

plot(family(kMember).traj(:,1), family(kMember).traj(:,3), ...
    'k-', 'LineWidth', 2.0);

for k = 1:nPhase
    for j = 1:nCircle
        X = unstableTraj{k,j}.x;
        plot(X(:,1), X(:,3), '-', ...
            'Color', [0.88 0.25 0.25], 'LineWidth', 0.45);

        X = stableTraj{k,j}.x;
        plot(X(:,1), X(:,3), '-', ...
            'Color', [0.25 0.72 0.40], 'LineWidth', 0.45);
    end
end

plot(-mu, 0, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
plot(1-mu, 0, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'k');
plot(L.L1(1), 0, 'rs', 'MarkerSize', 7, 'LineWidth', 1.4);

xlabel('$x~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('$z~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
title('Halo Manifold Tubes: $x$-$z$ Projection', ...
    'Interpreter', 'latex', 'FontSize', 15, 'FontWeight', 'bold');

% ----------------------------------------------------------
% Plot 3: x-y projection
% ----------------------------------------------------------
figure('Color','w');
hold on
grid on
box on
axis equal

plot(family(kMember).traj(:,1), family(kMember).traj(:,2), ...
    'k-', 'LineWidth', 2.0);

for k = 1:nPhase
    for j = 1:nCircle
        X = unstableTraj{k,j}.x;
        plot(X(:,1), X(:,2), '-', ...
            'Color', [0.88 0.25 0.25], 'LineWidth', 0.45);

        X = stableTraj{k,j}.x;
        plot(X(:,1), X(:,2), '-', ...
            'Color', [0.25 0.72 0.40], 'LineWidth', 0.45);
    end
end

plot(-mu, 0, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
plot(1-mu, 0, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'k');
plot(L.L1(1), 0, 'rs', 'MarkerSize', 7, 'LineWidth', 1.4);

xlabel('$x~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('$y~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
title('Halo Manifold Tubes: $x$-$y$ Projection', ...
    'Interpreter', 'latex', 'FontSize', 15, 'FontWeight', 'bold');

fprintf('\nDone.\n');