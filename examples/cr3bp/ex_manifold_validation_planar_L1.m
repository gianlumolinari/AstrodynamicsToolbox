clc;
clear;
close all;

% Validation example using a known corrected planar L1 Lyapunov orbit.

mu = 0.012150585609624;

state0 = [ ...
    0.8233908198365149;
    0;
    0;
    0;
    0.1263427298388180;
    0];

T = 2.7430007981241529;

orbit = astro.periodic.buildOrbitFromState(state0, T, mu, struct('verbose', true));

eps0 = 1e-6;
tfManifold = 3.0;

atlas = astro.manifolds.generateManifoldAtlas(orbit, mu, 'unstable', eps0, tfManifold, ...
    struct('numPhaseSamples', 40, ...
           'includeBothSigns', true, ...
           'normalizeMode', 'position', ...
           'RelTol', 1e-12, ...
           'AbsTol', 1e-12));

nB = numel(atlas.branches);
reports = struct( ...
    'type', {}, ...
    'C0', {}, ...
    'CmaxError', {}, ...
    'initialDistanceToReference', {}, ...
    'distanceAfterFirstSamples', {}, ...
    'appearsToDepart', {});
maxCerr = 0;
for k = 1:nB
    reports(k) = astro.manifolds.validateManifoldBranch(atlas.branches(k), orbit, mu);
    maxCerr = max(maxCerr, reports(k).CmaxError);
end

fprintf('\nValidation summary\n');
fprintf('------------------\n');
fprintf('Number of branches: %d\n', nB);
fprintf('Maximum Jacobi drift: %.3e\n', maxCerr);
fprintf('Orbit closure error: %.3e\n', orbit.closureError);

nDepart = 0;
for k = 1:nB
    if reports(k).appearsToDepart
        nDepart = nDepart + 1;
    end
end
fprintf('Branches initially departing from reference orbit: %d / %d\n', nDepart, nB);

lambda = atlas.eigData.eigvals;
fprintf('Largest |lambda|  : %.6e\n', max(abs(lambda)));
fprintf('Smallest |lambda| : %.6e\n', min(abs(lambda)));

fig = figure;
ax = axes(fig);
hold(ax,'on');
grid(ax,'on');

hU = astro.plot.plotManifoldFamily2D(atlas, struct( ...
    'ax', ax, ...
    'Color', 'c', ...
    'LineWidth', 0.8, ...
    'ShowOrbit', false, ...
    'ShowPrimaries', false));

hOrbit = plot(ax, orbit.x(:,1), orbit.x(:,2), 'k', 'LineWidth', 2.0);

hPrimary1 = plot(ax, -mu, 0, 'ko', 'MarkerFaceColor','k', 'MarkerSize', 6);
hPrimary2 = plot(ax, 1-mu, 0, 'o', ...
    'Color','r', ...
    'MarkerFaceColor',[0.7 0.7 0.7], ...
    'MarkerSize', 8);

L = astro.cr3bp.lagrangePoints(mu);
hL1 = plot(ax, L.L1(1), L.L1(2), 'ko', 'MarkerFaceColor','y', 'MarkerSize', 7);
hL2 = plot(ax, L.L2(1), L.L2(2), 'ko', 'MarkerFaceColor','c', 'MarkerSize', 7);

text(ax, L.L1(1)+0.004, L.L1(2)+0.006, '$L_1$', ...
    'Interpreter','latex', 'FontSize', 14);
text(ax, L.L2(1)+0.004, L.L2(2)+0.006, '$L_2$', ...
    'Interpreter','latex', 'FontSize', 14);

axis(ax,'equal');
xlabel(ax,'x');
ylabel(ax,'y');
title(ax,'Validation case: planar L1 unstable manifold');

legend(ax, ...
    [hU.branch(1), hOrbit, hPrimary1, hPrimary2, hL1, hL2], ...
    {'Unstable manifold', 'Periodic orbit', 'Earth', 'Moon', 'L_1', 'L_2'}, ...
    'Location','best');

xf = orbit.x(end,:).';
dx = xf - orbit.state0;

fprintf('State closure mismatch:\n');
fprintf('  dx  = %.3e\n', dx(1));
fprintf('  dy  = %.3e\n', dx(2));
fprintf('  dz  = %.3e\n', dx(3));
fprintf('  dvx = %.3e\n', dx(4));
fprintf('  dvy = %.3e\n', dx(5));
fprintf('  dvz = %.3e\n', dx(6));