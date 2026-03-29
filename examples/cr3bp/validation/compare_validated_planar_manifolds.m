clc;
clear;
close all;

% COMPARE_VALIDATED_PLANAR_MANIFOLDS
% Load saved validated planar Lyapunov members and compare their manifolds.

mu = 0.012150585609624;

orbitNames = { ...
    'planar_lyapunov_L1', ...
    'planar_lyapunov_L1_small', ...
    'planar_lyapunov_L1_medium', ...
    'planar_lyapunov_L1_large'};

eps0 = 2e-5;
tfStable = 3.0;
tfUnstable = 3.0;
numPhaseSamples = 30;

nOrb = numel(orbitNames);
manifolds = cell(nOrb,1);
orbits = cell(nOrb,1);

fprintf('\nComparing validated planar manifold cases\n');
fprintf('----------------------------------------\n');

for k = 1:nOrb
    S = astro.periodic.loadValidatedOrbit(orbitNames{k});

    orbit = astro.periodic.buildOrbitStruct(S.state0, S.T, S.mu, ...
        struct('family', S.family, ...
               'libr', S.libr, ...
               'system', 'earth-moon', ...
               'source', S.source, ...
               'dimension', '2D'), ...
        struct('verbose', true));

    man = astro.manifolds.generatePeriodicOrbitManifolds(orbit, S.mu, ...
        struct('eps0', eps0, ...
               'tfStable', tfStable, ...
               'tfUnstable', tfUnstable, ...
               'numPhaseSamples', numPhaseSamples, ...
               'includeBothSigns', true, ...
               'normalizeMode', 'position', ...
               'validateOrbit', true, ...
               'verbose', true));

    manifolds{k} = man;
    orbits{k} = orbit;

    fprintf('\nCase: %s\n', orbitNames{k});
    fprintf('  closure error            = %.3e\n', orbit.closureError);
    fprintf('  max stable Jacobi drift   = %.3e\n', man.maxJacobiDriftStable);
    fprintf('  max unstable Jacobi drift = %.3e\n', man.maxJacobiDriftUnstable);
end

%% ------------------------------------------------------------------------
% Orbit comparison
% -------------------------------------------------------------------------
fig1 = figure;
ax1 = axes(fig1);
hold(ax1,'on');
grid(ax1,'on');
box(ax1,'on');

for k = 1:nOrb
    X = orbits{k}.x;
    plot(ax1, X(:,1), X(:,2), 'LineWidth', 2.0);
end

L = astro.cr3bp.lagrangePoints(mu);
plot(ax1, -mu, 0, 'ko', 'MarkerFaceColor','k', 'MarkerSize', 6);
plot(ax1, 1-mu, 0, 'o', ...
    'Color',[0.4 0.4 0.4], ...
    'MarkerFaceColor',[0.7 0.7 0.7], ...
    'MarkerSize', 8);
plot(ax1, L.L1(1), L.L1(2), 'yo', 'MarkerFaceColor','y', 'MarkerSize', 7);

axis(ax1,'equal');
xlabel(ax1,'x');
ylabel(ax1,'y');
title(ax1,'Validated planar Lyapunov orbit database comparison');
legend(ax1, [orbitNames, {'Earth','Moon','L_1'}], 'Interpreter','none', 'Location','best');

%% ------------------------------------------------------------------------
% One manifold figure per saved member
% -------------------------------------------------------------------------
for k = 1:nOrb
    man = manifolds{k};
    orbit = orbits{k};

    fig = figure;
    ax = axes(fig);
    hold(ax,'on');
    grid(ax,'on');
    box(ax,'on');

    hU = astro.plot.plotManifoldFamily2D(man.unstable, struct( ...
        'ax', ax, ...
        'Color', [1 0 0], ...
        'LineWidth', 0.8, ...
        'ShowOrbit', false, ...
        'ShowPrimaries', false));

    hS = astro.plot.plotManifoldFamily2D(man.stable, struct( ...
        'ax', ax, ...
        'Color', [0 0 1], ...
        'LineWidth', 0.8, ...
        'ShowOrbit', false, ...
        'ShowPrimaries', false));

    hOrbit = plot(ax, orbit.x(:,1), orbit.x(:,2), 'k', 'LineWidth', 2.0);

    hEarth = plot(ax, -mu, 0, 'ko', 'MarkerFaceColor','k', 'MarkerSize', 6);
    hMoon  = plot(ax, 1-mu, 0, 'o', ...
        'Color',[0.4 0.4 0.4], ...
        'MarkerFaceColor',[0.7 0.7 0.7], ...
        'MarkerSize', 8);
    hL1 = plot(ax, L.L1(1), L.L1(2), 'yo', 'MarkerFaceColor','y', 'MarkerSize', 7);

    axis(ax,'equal');
    xlabel(ax,'x');
    ylabel(ax,'y');
    title(ax, sprintf('Stable/unstable manifolds: %s', orbitNames{k}), ...
        'Interpreter','none');

    legend(ax, ...
        [hU.branch(1), hS.branch(1), hOrbit, hEarth, hMoon, hL1], ...
        {'Unstable manifold', 'Stable manifold', 'Periodic orbit', ...
         'Earth', 'Moon', 'L_1'}, ...
        'Location', 'best');
end