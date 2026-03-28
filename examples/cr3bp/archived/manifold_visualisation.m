clc;
clear;
close all;

% =========================================================================%
% This script:
%   1) queries the JPL periodic orbit API for Earth-Moon planar L1 Lyapunov
%      seeds in a small-period range
%   2) tests several candidate rows
%   3) locally corrects each candidate with your planar Lyapunov corrector
%   4) selects the corrected orbit closest to the local L1 region
%   5) generates stable and unstable manifolds
%   6) plots them in Moon-centered km
%
% REQUIREMENTS
%   - astro.cr3bp.queryJPLPeriodicOrbits
%   - astro.cr3bp.differentialCorrectionPlanarLyapunov
%   - astro.cr3bp.lagrangePoints
%   - astro.periodic.buildOrbitFromState
%   - astro.manifolds.generateManifoldAtlas
%   - astro.manifolds.sortTubeStrands
%   - astro.manifolds.filterBranchesBySign   (optional but recommended)
% =========================================================================

%% User settings
mu = 0.012150585609624;      % Earth-Moon
libr = 1;                    % target L1
periodMin = 2.6;             % TU
periodMax = 2.8;             % TU
nCandidateRows = 20;         % test first N rows returned by JPL
maxIterCorr = 30;
tolCorr = 1e-12;

% Manifold settings
eps0 = 1e-5;
tfUnstable = 5.0;
tfStable   = 5.0;
numSamples = 50;
plotOnlyOneSignEach = false; 

figureName = 'Planar L1 Lyapunov manifolds';

%% Lagrange points in normalized barycentric frame
L = astro.cr3bp.lagrangePoints(mu);
xL1 = L.L1(1);
xMoon = 1 - mu;

fprintf('Computed normalized L1 location: %.16f\n', xL1);

%% ------------------------------------------------------------------------
% Query JPL planar L1 Lyapunov seeds in a small-period range
% -------------------------------------------------------------------------
data = astro.cr3bp.queryJPLPeriodicOrbits( ...
    'sys', 'earth-moon', ...
    'family', 'lyapunov', ...
    'libr', libr, ...
    'periodmin', periodMin, ...
    'periodmax', periodMax, ...
    'periodunits', 'TU');

fields = string(data.fields);
rows = data.data;

nRows = min(nCandidateRows, numel(rows));
fprintf('Testing %d JPL candidate rows in period range [%.2f, %.2f] TU\n', ...
    nRows, periodMin, periodMax);

LU = localExtractLengthUnitKm(data);   % km
fprintf('Using length unit LU = %.3f km\n', LU);

%% ------------------------------------------------------------------------
% Scan candidate rows, correct locally, keep the orbit closest to local L1
% -------------------------------------------------------------------------
best = struct();
best.found = false;
best.metric = inf;

for i = 1:nRows
    try
        row = localNormalizeRow(rows{i});

        seed.x0 = localGetFieldValue(row, fields, 'x');
        seed.y0 = localGetFieldValue(row, fields, 'y');
        seed.z0 = localGetFieldValue(row, fields, 'z');
        seed.vx0 = localGetFieldValue(row, fields, 'vx');
        seed.vy0 = localGetFieldValue(row, fields, 'vy');
        seed.vz0 = localGetFieldValue(row, fields, 'vz');
        seed.C  = localGetFieldValue(row, fields, 'jacobi');
        seed.T  = localGetFieldValue(row, fields, 'period');

        % Only use rows already close to planar symmetry form
        planarResidual = norm([seed.y0, seed.z0, seed.vx0, seed.vz0]);
        if planarResidual > 1e-6
            fprintf('  Row %3d skipped: not close to planar symmetry form (%.3e)\n', ...
                i, planarResidual);
            continue
        end

        fprintf('\nCandidate row %d\n', i);
        fprintf('  seed x0   = %.16f\n', seed.x0);
        fprintf('  seed vy0  = %.16f\n', seed.vy0);
        fprintf('  seed T    = %.16f\n', seed.T);
        fprintf('  seed C    = %.16f\n', seed.C);

        corr = astro.cr3bp.differentialCorrectionPlanarLyapunov( ...
            seed.x0, seed.vy0, mu, maxIterCorr, tolCorr);

        if ~corr.converged
            fprintf('  -> correction failed to converge\n');
            continue
        end

        state0_i = corr.state0;
        T_i = corr.period;
        orbit_i = astro.periodic.buildOrbitFromState(state0_i, T_i, mu, struct('verbose', false));

        dxL1 = abs(state0_i(1) - xL1);

        fprintf('  corrected x0 = %.16f\n', state0_i(1));
        fprintf('  corrected vy0= %.16f\n', state0_i(5));
        fprintf('  corrected T  = %.16f\n', T_i);
        fprintf('  closure err  = %.3e\n', orbit_i.closureError);
        fprintf('  |x0 - xL1|   = %.3e\n', dxL1);

        % Prefer:
        %   - small closure error
        %   - orbit close to local L1 region
        metric = dxL1 + 1e3 * orbit_i.closureError;

        if orbit_i.closureError < 1e-8 && metric < best.metric
            best.found = true;
            best.metric = metric;
            best.rowIndex = i;
            best.seed = seed;
            best.corr = corr;
            best.state0 = state0_i;
            best.T = T_i;
            best.orbit = orbit_i;
            best.dxL1 = dxL1;
        end

    catch ME
        fprintf('  Row %d failed: %s\n', i, ME.message);
    end
end

if ~best.found
    error('No suitable corrected planar L1 orbit found in the tested JPL rows.');
end

%% Selected orbit
state0 = best.state0;
T = best.T;
orbit = best.orbit;

fprintf('\n============================================================\n');
fprintf('Selected candidate row: %d\n', best.rowIndex);
fprintf('Selected corrected orbit:\n');
fprintf('  x0         = %.16f\n', state0(1));
fprintf('  vy0        = %.16f\n', state0(5));
fprintf('  T          = %.16f\n', T);
fprintf('  closure    = %.3e\n', orbit.closureError);
fprintf('  |x0 - xL1| = %.3e\n', best.dxL1);
fprintf('============================================================\n');

%% Save validated orbit
if ~exist('data/validated', 'dir')
    mkdir('data/validated');
end
family = 'lyapunov';
source = sprintf('JPL periodic orbit API row %d corrected locally', best.rowIndex);
closureError = orbit.closureError;
save('data/validated/planarLyapunov_L1_Candidate.mat', ...
    'state0', 'T', 'mu', 'family', 'libr', 'source', 'closureError');

%% ------------------------------------------------------------------------
% Generate manifolds
% -------------------------------------------------------------------------
atlasU = astro.manifolds.generateManifoldAtlas(orbit, mu, 'unstable', eps0, tfUnstable, ...
    struct('numPhaseSamples', numSamples, ...
           'includeBothSigns', true, ...
           'normalizeMode', 'position', ...
           'RelTol', 1e-12, ...
           'AbsTol', 1e-12));

atlasS = astro.manifolds.generateManifoldAtlas(orbit, mu, 'stable', eps0, tfStable, ...
    struct('numPhaseSamples', numSamples, ...
           'includeBothSigns', true, ...
           'normalizeMode', 'position', ...
           'RelTol', 1e-12, ...
           'AbsTol', 1e-12));

atlasU = astro.manifolds.sortTubeStrands(atlasU, 'signphase');
atlasS = astro.manifolds.sortTubeStrands(atlasS, 'signphase');

if plotOnlyOneSignEach
    atlasU = astro.manifolds.filterBranchesBySign(atlasU, +1);
    atlasS = astro.manifolds.filterBranchesBySign(atlasS, -1);
end

%% Diagnostics
cerrU = max(arrayfun(@(b) b.CerrMax, atlasU.branches));
cerrS = max(arrayfun(@(b) b.CerrMax, atlasS.branches));

fprintf('\nManifold summary\n');
fprintf('----------------\n');
fprintf('Unstable branches: %d\n', numel(atlasU.branches));
fprintf('Stable branches  : %d\n', numel(atlasS.branches));
fprintf('Max Jacobi drift (unstable): %.3e\n', cerrU);
fprintf('Max Jacobi drift (stable)  : %.3e\n', cerrS);

%% ------------------------------------------------------------------------
% Convert to Moon-centered km for plotting
% -------------------------------------------------------------------------
[xOrbit_km, yOrbit_km] = localToMoonCenteredKm(orbit.x(:,1), orbit.x(:,2), xMoon, LU);
[xL1_km, yL1_km]       = localToMoonCenteredKm(L.L1(1), L.L1(2), xMoon, LU);
[xL2_km, yL2_km]       = localToMoonCenteredKm(L.L2(1), L.L2(2), xMoon, LU);

%% ------------------------------------------------------------------------
% Plot
% -------------------------------------------------------------------------
fig = figure('Color', 'w');
ax = axes(fig);
hold(ax, 'on');
grid(ax, 'on');
box(ax, 'on');

% Stable manifolds
for k = 1:numel(atlasS.branches)
    X = atlasS.branches(k).x;
    [xkm, ykm] = localToMoonCenteredKm(X(:,1), X(:,2), xMoon, LU);
    plot(ax, xkm, ykm, 'g', 'LineWidth', 1);
end

% Unstable manifolds
for k = 1:numel(atlasU.branches)
    X = atlasU.branches(k).x;
    [xkm, ykm] = localToMoonCenteredKm(X(:,1), X(:,2), xMoon, LU);
    plot(ax, xkm, ykm, 'r', 'LineWidth', 1);
end

% Periodic orbit in black
plot(ax, xOrbit_km, yOrbit_km, 'k', 'LineWidth', 2.0);

% Moon at origin
plot(ax, 0, 0, 'o', ...
    'MarkerFaceColor', [0.6 0.6 0.6], ...
    'MarkerEdgeColor', [0.4 0.4 0.4], ...
    'MarkerSize', 10);

% L1 point
plot(ax, xL1_km, yL1_km, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);

xlabel(ax, '$x\ (\times 10^4\ \mathrm{km})$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel(ax, '$y\ (\times 10^4\ \mathrm{km})$', 'Interpreter', 'latex', 'FontSize', 18);
title(ax, figureName, 'Interpreter', 'latex');

% Scale to x10^4 km visually
xt = get(ax, 'XTick');
yt = get(ax, 'YTick');
set(ax, 'XTickLabel', compose('%.0f', xt/1e4));
set(ax, 'YTickLabel', compose('%.0f', yt/1e4));

axis(ax, 'equal');

% xlim(ax, [-14e4, 2.5e4]);
% ylim(ax, [-7e4, 7e4]);

set(ax, 'FontSize', 16);

% Labels
text(ax, xL1_km + 0.1e4, yL1_km + 2.0e4, '$L_1$', ...
    'Interpreter', 'latex', 'FontSize', 18);
text(ax, 0.2e4, 2.5e4, 'Moon', 'Interpreter', 'latex', 'FontSize', 18);



%% ------------------------------------------------------------------------
% Helper functions
% -------------------------------------------------------------------------
function row = localNormalizeRow(rawRow)
    row = rawRow;
    while iscell(row) && numel(row) == 1
        row = row{1};
    end
    row = row(:).';
end

function val = localGetFieldValue(row, fields, targetName)
    idx = find(strcmpi(fields, targetName), 1);
    raw = row{idx};

    while iscell(raw) && numel(raw) == 1
        raw = raw{1};
    end

    if isnumeric(raw)
        val = double(raw);
    elseif isstring(raw) || ischar(raw)
        val = str2double(raw);
    else
        error('Unsupported field type for "%s".', targetName);
    end
end

function LU = localExtractLengthUnitKm(data)
    LU = 384400.0;
    if isstruct(data)
        if isfield(data, 'system') && isstruct(data.system)
            if isfield(data.system, 'lunit')
                raw = data.system.lunit;
                if isnumeric(raw)
                    LU = double(raw);
                    return
                elseif ischar(raw) || isstring(raw)
                    tmp = str2double(raw);
                    if ~isnan(tmp), LU = tmp; return; end
                end
            end
        end
    end
end

function [xkm, ykm] = localToMoonCenteredKm(x, y, xMoon, LU)
    xkm = (x - xMoon) * LU;
    ykm = y * LU;
end