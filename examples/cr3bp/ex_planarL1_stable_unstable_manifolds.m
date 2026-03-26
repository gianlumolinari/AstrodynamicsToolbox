clc;
clear;
close all;


mu = 0.012150585609624;

state0 = [0.8233908198365149;
    0;
    0;
    0;
    0.1263427298388180;
    0];

T = 2.7430007981241529;

orbit = astro.periodic.buildOrbitFromState(state0, T, mu, struct('verbose',true));

% -------------------------------------------------------------------------
% Manifold settings
% -------------------------------------------------------------------------
eps0 = 1e-5;
tfUnstable = 5.0;
tfStable   = 5.0;

numSamples = 30;

atlasU = astro.manifolds.generateManifoldAtlas(orbit, mu, 'unstable', eps0, tfUnstable, ...
    struct('numPhaseSamples',numSamples, ...
           'includeBothSigns',true, ...
           'normalizeMode','position', ...
           'RelTol',1e-12, ...
           'AbsTol',1e-12));

atlasS = astro.manifolds.generateManifoldAtlas(orbit, mu, 'stable', eps0, tfStable, ...
    struct('numPhaseSamples',numSamples, ...
           'includeBothSigns',true, ...
           'normalizeMode','position', ...
           'RelTol',1e-12, ...
           'AbsTol',1e-12));

atlasU = astro.manifolds.sortTubeStrands(atlasU, 'signphase');
atlasS = astro.manifolds.sortTubeStrands(atlasS, 'signphase');


fig = figure;
ax = axes(fig);
hold(ax,'on');
grid(ax,'on');

hU = astro.plot.plotManifoldFamily2D(atlasU, struct( ...
    'ax', ax, ...
    'Color', 'r', ...
    'LineWidth', 0.8, ...
    'ShowOrbit', false, ...
    'ShowPrimaries', false));

hS = astro.plot.plotManifoldFamily2D(atlasS, struct( ...
    'ax', ax, ...
    'Color', 'g', ...
    'LineWidth', 0.8, ...
    'ShowOrbit', false, ...
    'ShowPrimaries', false));

% Plot periodic orbit last
hOrbit = plot(ax, orbit.x(:,1), orbit.x(:,2), 'k', 'LineWidth', 2.0);

% Primaries
hPrimary1 = plot(ax, -mu, 0, 'ko', 'MarkerFaceColor','k', 'MarkerSize', 6);
hPrimary2 = plot(ax, 1-mu, 0, 'o', ...
    'Color',[0.4 0.4 0.4], ...
    'MarkerFaceColor',[0.7 0.7 0.7], ...
    'MarkerSize', 8);

axis(ax,'equal');
xlabel(ax,'x');
ylabel(ax,'y');
title(ax,'Planar L1 stable and unstable manifolds');

% L points
L = astro.cr3bp.lagrangePoints(mu);

hL1 = plot(ax, L.L1(1), L.L1(2), 'ko', 'MarkerFaceColor','y', 'MarkerSize', 7);
hL2 = plot(ax, L.L2(1), L.L2(2), 'ko', 'MarkerFaceColor','c', 'MarkerSize', 7);

text(ax, L.L1(1)+0.004, L.L1(2)+0.006, '$L_1$', ...
    'Interpreter','latex', 'FontSize', 14);

text(ax, L.L2(1)+0.004, L.L2(2)+0.006, '$L_2$', ...
    'Interpreter','latex', 'FontSize', 14);

xlim(ax, [-1 1.5]);
ylim(ax, [-1 1]);

set(ax,'FontSize',14);
box(ax,'on');

legend(ax, ...
    [hU.branch(1), hS.branch(1), hOrbit, hPrimary1, hPrimary2, hL1, hL2], ...
    {'Unstable manifold', 'Stable manifold', 'Periodic orbit', ...
     'Primary 1', 'Primary 2', 'L_1', 'L_2'}, ...
    'Location','best');

% -------------------------------------------------------------------------
% Diagnostics
% -------------------------------------------------------------------------
cerrU = max(arrayfun(@(b) b.CerrMax, atlasU.branches));
cerrS = max(arrayfun(@(b) b.CerrMax, atlasS.branches));

fprintf('Max Jacobi drift, unstable atlas: %.3e\n', cerrU);
fprintf('Max Jacobi drift, stable   atlas: %.3e\n', cerrS);