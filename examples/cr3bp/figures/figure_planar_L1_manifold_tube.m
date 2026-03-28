clc;
clear;
close all;

% EX_PLANAR_LYAPUNOV_MANIFOLD_DRIVER
% One unified script to:
%   1) fetch a JPL seed
%   2) correct it locally
%   3) save a validated periodic orbit
%   4) validate manifold propagation
%   5) generate stable/unstable manifold figures

mu = 0.012150585609624;

%% ========================================================================
% USER OPTIONS
% ========================================================================
libr = 1;              % 1 for L1, 2 for L2
jplIndex = 134;          % row index returned by JPL API
saveValidatedOrbit = true;

% Validation atlas settings
epsVal = 1e-6;
tfVal  = 3.0;
numSamplesVal = 40;

% Figure-generation atlas settings
epsFig = 2e-5;
tfUnstable = 5.0;
tfStable   = 5.0;
numSamplesFig = 25;

%% ========================================================================
% 1) FETCH JPL SEED
% ========================================================================
seed = astro.cr3bp.getPlanarLyapunovSeed('jpl', ...
    'sys', 'earth-moon', ...
    'family', 'lyapunov', ...
    'libr', libr, ...
    'index', jplIndex);

fprintf('Using seed from: %s\n', seed.source);
fprintf('  x0_seed = %.16f\n', seed.x0);
fprintf('  vy0_seed = %.16f\n', seed.vy0);
fprintf('  T_seed = %.16f\n', seed.T);
fprintf('  C_seed = %.16f\n', seed.C);

%% ========================================================================
% 2) LOCAL DIFFERENTIAL CORRECTION
% ========================================================================
corr = astro.cr3bp.differentialCorrectionPlanarLyapunov( ...
    seed.x0, seed.vy0, mu, 30, 1e-12);

if ~corr.converged
    error('Differential correction did not converge from the supplied seed.');
end

state0 = corr.state0;
T = corr.period;

fprintf('\nCorrected orbit\n');
fprintf('---------------\n');
fprintf('  x0 = %.16f\n', state0(1));
fprintf('  vy0 = %.16f\n', state0(5));
fprintf('  T  = %.16f\n', T);

orbit = astro.periodic.buildOrbitFromState(state0, T, mu, struct('verbose', true));

%% ========================================================================
% 3) SAVE VALIDATED ORBIT
% ========================================================================
validatedFile = sprintf('data/validated/planarLyapunov_L%d_validated.mat', libr);

if saveValidatedOrbit && orbit.closureError < 1e-8
    if ~exist('data/validated','dir')
        mkdir('data/validated');
    end

    family = 'lyapunov';
    source = seed.source;
    closureError = orbit.closureError;
    seedIndex = jplIndex;

    save(validatedFile, ...
        'state0', 'T', 'mu', 'family', 'libr', 'source', ...
        'closureError', 'seedIndex');

    fprintf('Saved corrected orbit to %s\n', validatedFile);
else
    fprintf('Orbit not saved: closure error is still too large or saving disabled.\n');
end

%% ========================================================================
% 4) VALIDATION ATLAS (UNSTABLE ONLY)
% ========================================================================
atlasVal = astro.manifolds.generateManifoldAtlas(orbit, mu, 'unstable', epsVal, tfVal, ...
    struct('numPhaseSamples', numSamplesVal, ...
           'includeBothSigns', true, ...
           'normalizeMode', 'position', ...
           'RelTol', 1e-12, ...
           'AbsTol', 1e-12));

nB = numel(atlasVal.branches);

reports = struct( ...
    'type', {}, ...
    'C0', {}, ...
    'CmaxError', {}, ...
    'initialDistanceToReference', {}, ...
    'distanceAfterFirstSamples', {}, ...
    'appearsToDepart', {});

maxCerrVal = 0;
nDepart = 0;

for k = 1:nB
    reports(k,1) = astro.manifolds.validateManifoldBranch(atlasVal.branches(k), orbit, mu);
    maxCerrVal = max(maxCerrVal, reports(k).CmaxError);
    if reports(k).appearsToDepart
        nDepart = nDepart + 1;
    end
end

lambda = atlasVal.eigData.eigvals;

fprintf('\nValidation summary\n');
fprintf('------------------\n');
fprintf('Number of validation branches: %d\n', nB);
fprintf('Maximum Jacobi drift: %.3e\n', maxCerrVal);
fprintf('Orbit closure error: %.3e\n', orbit.closureError);
fprintf('Branches initially departing from reference orbit: %d / %d\n', nDepart, nB);
fprintf('Largest |lambda|  : %.6e\n', max(abs(lambda)));
fprintf('Smallest |lambda| : %.6e\n', min(abs(lambda)));

xf = orbit.x(end,:).';
dx = xf - orbit.state0;

fprintf('State closure mismatch:\n');
fprintf('  dx  = %.3e\n', dx(1));
fprintf('  dy  = %.3e\n', dx(2));
fprintf('  dz  = %.3e\n', dx(3));
fprintf('  dvx = %.3e\n', dx(4));
fprintf('  dvy = %.3e\n', dx(5));
fprintf('  dvz = %.3e\n', dx(6));

%% ========================================================================
% 5) GENERATION ATLASES (STABLE + UNSTABLE)
% ========================================================================
atlasU = astro.manifolds.generateManifoldAtlas(orbit, mu, 'unstable', epsFig, tfUnstable, ...
    struct('numPhaseSamples', numSamplesFig, ...
           'includeBothSigns', true, ...
           'normalizeMode', 'position', ...
           'RelTol', 1e-12, ...
           'AbsTol', 1e-12));

atlasS = astro.manifolds.generateManifoldAtlas(orbit, mu, 'stable', epsFig, tfStable, ...
    struct('numPhaseSamples', numSamplesFig, ...
           'includeBothSigns', true, ...
           'normalizeMode', 'position', ...
           'RelTol', 1e-12, ...
           'AbsTol', 1e-12));

atlasU = astro.manifolds.sortTubeStrands(atlasU, 'signphase');
atlasS = astro.manifolds.sortTubeStrands(atlasS, 'signphase');

cerrU = max(arrayfun(@(b) b.CerrMax, atlasU.branches));
cerrS = max(arrayfun(@(b) b.CerrMax, atlasS.branches));

fprintf('\nManifold generation summary\n');
fprintf('---------------------------\n');
fprintf('Unstable branches: %d\n', numel(atlasU.branches));
fprintf('Stable branches  : %d\n', numel(atlasS.branches));
fprintf('Max Jacobi drift (unstable): %.3e\n', cerrU);
fprintf('Max Jacobi drift (stable)  : %.3e\n', cerrS);

%% ========================================================================
% 6) COMMON PLOT ELEMENTS
% ========================================================================
L = astro.cr3bp.lagrangePoints(mu);

plotLabelL1 = sprintf('$L_%d$', 1);
plotLabelL2 = sprintf('$L_%d$', 2);
plotTitleBase = sprintf('Planar L%d Lyapunov manifolds', libr);

%% ========================================================================
% 7) VALIDATION FIGURE
% ========================================================================
figVal = figure;
axVal = axes(figVal);
hold(axVal,'on');
grid(axVal,'on');

hVal = astro.plot.plotManifoldFamily2D(atlasVal, struct( ...
    'ax', axVal, ...
    'Color', 'c', ...
    'LineWidth', 0.8, ...
    'ShowOrbit', false, ...
    'ShowPrimaries', false));

hOrbitVal = plot(axVal, orbit.x(:,1), orbit.x(:,2), 'k', 'LineWidth', 2.0);
hEarthVal = plot(axVal, -mu, 0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
hMoonVal  = plot(axVal, 1-mu, 0, 'o', 'Color', [0.4 0.4 0.4], ...
    'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerSize', 8);

hL1Val = plot(axVal, L.L1(1), L.L1(2), 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 7);
hL2Val = plot(axVal, L.L2(1), L.L2(2), 'ko', 'MarkerFaceColor', 'c', 'MarkerSize', 7);

text(axVal, L.L1(1)+0.004, L.L1(2)+0.006, plotLabelL1, ...
    'Interpreter', 'latex', 'FontSize', 14);
text(axVal, L.L2(1)+0.004, L.L2(2)+0.006, plotLabelL2, ...
    'Interpreter', 'latex', 'FontSize', 14);

axis(axVal,'equal');
xlabel(axVal,'x');
ylabel(axVal,'y');
title(axVal, sprintf('Validation: planar L%d unstable manifold', libr));

legend(axVal, ...
    [hVal.branch(1), hOrbitVal, hEarthVal, hMoonVal, hL1Val, hL2Val], ...
    {'Unstable manifold', 'Periodic orbit', 'Earth', 'Moon', 'L_1', 'L_2'}, ...
    'Location', 'best');

set(axVal,'FontSize',14);
box(axVal,'on');

%% ========================================================================
% 8) COMBINED STABLE + UNSTABLE FIGURE
% ========================================================================
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

hOrbit = plot(ax, orbit.x(:,1), orbit.x(:,2), 'k', 'LineWidth', 2.0);
hEarth = plot(ax, -mu, 0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
hMoon  = plot(ax, 1-mu, 0, 'o', 'Color', [0.4 0.4 0.4], ...
    'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerSize', 8);

hL1 = plot(ax, L.L1(1), L.L1(2), 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 7);
hL2 = plot(ax, L.L2(1), L.L2(2), 'ko', 'MarkerFaceColor', 'c', 'MarkerSize', 7);

text(ax, L.L1(1)+0.004, L.L1(2)+0.006, plotLabelL1, ...
    'Interpreter', 'latex', 'FontSize', 14);
text(ax, L.L2(1)+0.004, L.L2(2)+0.006, plotLabelL2, ...
    'Interpreter', 'latex', 'FontSize', 14);

axis(ax,'equal');
xlabel(ax,'x');
ylabel(ax,'y');
title(ax, plotTitleBase);

set(ax,'FontSize',14);
box(ax,'on');

legend(ax, ...
    [hU.branch(1), hS.branch(1), hOrbit, hEarth, hMoon, hL1, hL2], ...
    {'Unstable manifold', 'Stable manifold', 'Periodic orbit', ...
     'Earth', 'Moon', 'L_1', 'L_2'}, ...
    'Location', 'best');

%% ========================================================================
% 9) OPTIONAL UNSTABLE-ONLY FIGURE
% ========================================================================
fig2 = figure;
ax2 = axes(fig2);
hold(ax2,'on');
grid(ax2,'on');

hU2 = astro.plot.plotManifoldFamily2D(atlasU, struct( ...
    'ax', ax2, ...
    'Color', 'c', ...
    'LineWidth', 0.8, ...
    'ShowOrbit', false, ...
    'ShowPrimaries', false));

hOrbit2 = plot(ax2, orbit.x(:,1), orbit.x(:,2), 'k', 'LineWidth', 2.0);
plot(ax2, -mu, 0, 'ko', 'MarkerFaceColor','k', 'MarkerSize', 6);
plot(ax2, 1-mu, 0, 'o', 'Color',[0.4 0.4 0.4], ...
    'MarkerFaceColor',[0.7 0.7 0.7], 'MarkerSize', 8);
plot(ax2, L.L1(1), L.L1(2), 'ko', 'MarkerFaceColor','y', 'MarkerSize', 7);
plot(ax2, L.L2(1), L.L2(2), 'ko', 'MarkerFaceColor','c', 'MarkerSize', 7);

text(ax2, L.L1(1)+0.004, L.L1(2)+0.006, plotLabelL1, ...
    'Interpreter', 'latex', 'FontSize', 14);
text(ax2, L.L2(1)+0.004, L.L2(2)+0.006, plotLabelL2, ...
    'Interpreter', 'latex', 'FontSize', 14);

axis(ax2,'equal');
xlabel(ax2,'x');
ylabel(ax2,'y');
title(ax2, sprintf('Planar L%d unstable manifold', libr));

legend(ax2, [hU2.branch(1), hOrbit2], ...
    {'Unstable manifold', 'Periodic orbit'}, ...
    'Location', 'best');