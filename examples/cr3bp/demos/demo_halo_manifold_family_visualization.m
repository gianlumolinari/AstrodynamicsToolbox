clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/AstroToolbox')
startup

% ==========================================================
% DEMO_HALO_MANIFOLD_FAMILY_VISUALIZATION
% Earth-Moon L1 northern halo family
%
% This version is adapted to the current toolbox structure:
%   - halo family from JPL seeds + pseudo-arclength continuation
%   - orbit struct built with astro.periodic.buildOrbitStruct
%   - eigenstructure from astro.manifolds.computeEigenstructure
%   - manifold seeds from astro.manifolds.makeManifoldSeeds
%   - branch propagation done directly with astro.propagators.propagate
%
% NOTE:
% This visualizes stable/unstable manifold branches seeded along the orbit.
% It is not the old ring-based "tube surface" construction.
% ==========================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('Halo Manifold Family Visualization Demo\n');
fprintf('============================================================\n');

%% ----------------------------------------------------------
% 1) Build halo family and pick one member
% -----------------------------------------------------------
data = astro.cr3bp.queryJPLPeriodicOrbits( ...
    'sys', 'earth-moon', ...
    'family', 'halo', ...
    'libr', 1, ...
    'branch', 'N');

seed1 = astro.cr3bp.parseJPLPeriodicOrbit(data, 1);
seed2 = astro.cr3bp.parseJPLPeriodicOrbit(data, 2);

mu = seed1.mu;

family = astro.cr3bp.continueHaloPseudoArc(seed1, seed2, 20, 5e-4, mu);

kMember = 18;
if kMember > numel(family)
    error('Requested halo family member exceeds generated family size.');
end

state0 = family(kMember).state0(:);
T      = family(kMember).period;
C      = family(kMember).C;

fprintf('\nSelected halo family member %d\n', kMember);
fprintf('  Period = %.12f TU\n', T);
fprintf('  Jacobi = %.12f\n', C);

%% ----------------------------------------------------------
% 2) Build orbit struct and eigenspace data
% -----------------------------------------------------------
orbit = astro.periodic.buildOrbitStruct(state0, T, mu, ...
    struct('family','halo', ...
           'libr',1, ...
           'system','earth-moon', ...
           'source','demo_halo_manifold_family_visualization', ...
           'dimension','3D'), ...
    struct('verbose', true));

eigData = astro.manifolds.computeEigenstructure(orbit, mu, struct('verbose', true));

%% ----------------------------------------------------------
% 3) Seed settings
% -----------------------------------------------------------
nPhase = 15;
epsManifold = 1e-5;

seedOpts = struct();
seedOpts.numPhaseSamples = nPhase;
seedOpts.includeBothSigns = true;

unstableSeeds = astro.manifolds.makeManifoldSeeds( ...
    orbit, eigData, 'unstable', epsManifold, seedOpts);

stableSeeds = astro.manifolds.makeManifoldSeeds( ...
    orbit, eigData, 'stable', epsManifold, seedOpts);

fprintf('\nGenerated manifold seeds:\n');
fprintf('  Unstable seeds : %d\n', numel(unstableSeeds));
fprintf('  Stable seeds   : %d\n', numel(stableSeeds));

%% ----------------------------------------------------------
% 4) Propagation settings
% -----------------------------------------------------------
tUnstable = 6.0;
tStable   = -6.0;

odeOpts = struct();
odeOpts.RelTol = 1e-12;
odeOpts.AbsTol = 1e-12;
odeOpts.Solver = 'ode113';

unstableTraj = cell(numel(unstableSeeds),1);
stableTraj   = cell(numel(stableSeeds),1);

fprintf('\nPropagating unstable manifold branches...\n');
for k = 1:numel(unstableSeeds)
    unstableTraj{k} = astro.propagators.propagate( ...
        @(t,x) astro.cr3bp.eomCR3BP(t,x,mu), ...
        [0, tUnstable], unstableSeeds(k).x0, odeOpts);
end

fprintf('Propagating stable manifold branches...\n');
for k = 1:numel(stableSeeds)
    stableTraj{k} = astro.propagators.propagate( ...
        @(t,x) astro.cr3bp.eomCR3BP(t,x,mu), ...
        [0, tStable], stableSeeds(k).x0, odeOpts);
end

fprintf('Max unstable final distance from orbit seed set: %.3e\n', ...
    max(cellfun(@(s) norm(s.x(end,:)' - s.x(1,:)'), unstableTraj)));

fprintf('Max stable final distance from orbit seed set: %.3e\n', ...
    max(cellfun(@(s) norm(s.x(end,:)' - s.x(1,:)'), stableTraj)));

%% ----------------------------------------------------------
% 5) Lagrange points / primaries
% -----------------------------------------------------------
L = astro.cr3bp.lagrangePoints(mu);

[xL1, yL1, zL1, xL2, yL2, zL2] = localParseLagrangeOutput(L);

%% ----------------------------------------------------------
% 6) Plot 1: 3D halo and manifold families
% -----------------------------------------------------------
fig1 = figure;
ax1 = axes(fig1);
hold(ax1,'on');
grid(ax1,'on');
box(ax1,'on');
axis(ax1,'equal');

hOrbit3D = plot3(ax1, orbit.x(:,1), orbit.x(:,2), orbit.x(:,3), ...
    'y-', 'LineWidth', 2.0);

hU3D = gobjects(numel(unstableTraj),1);
for k = 1:numel(unstableTraj)
    X = unstableTraj{k}.x;
    hU3D(k) = plot3(ax1, X(:,1), X(:,2), X(:,3), '-', ...
        'Color', [0.88 0.25 0.25], 'LineWidth', 1);
end

hS3D = gobjects(numel(stableTraj),1);
for k = 1:numel(stableTraj)
    X = stableTraj{k}.x;
    hS3D(k) = plot3(ax1, X(:,1), X(:,2), X(:,3), '-', ...
        'Color', [0.25 0.72 0.40], 'LineWidth', 1);
end

hPrimary = plot3(ax1, -mu, 0, 0, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
hSecondary = plot3(ax1, 1-mu, 0, 0, 'co', 'MarkerSize', 7, 'MarkerFaceColor', 'k');
hL1 = plot3(ax1, xL1, yL1, zL1, 'md', 'MarkerSize', 7, 'LineWidth', 1.4);
hL2 = plot3(ax1, xL2, yL2, zL2, 'md', 'MarkerSize', 7, 'LineWidth', 1.4);

xlabel(ax1, '$x~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
ylabel(ax1, '$y~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
zlabel(ax1, '$z~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
title(ax1, sprintf('Halo Orbit with Stable/Unstable Manifold Families $(C = %.6f)$', C), ...
    'Interpreter', 'latex', 'FontSize', 15, 'FontWeight', 'bold');
legend(ax1, [hOrbit3D, hU3D(1), hS3D(1), hPrimary, hSecondary, hL1, hL2], ...
    {'Halo orbit', 'Unstable family', 'Stable family', 'Primary', 'Secondary', '$L_1$', '$L_2$'}, ...
    'Interpreter', 'latex', 'Location', 'best');
view(3);

fprintf('\nDone.\n');

%% =========================================================
% Local helper
% ==========================================================
function [xL1, yL1, zL1, xL2, yL2, zL2] = localParseLagrangeOutput(L)
    if isstruct(L)
        xL1 = L.L1(1); yL1 = L.L1(2); zL1 = L.L1(3);
        xL2 = L.L2(1); yL2 = L.L2(2); zL2 = L.L2(3);
    elseif isnumeric(L) && size(L,2) == 3 && size(L,1) >= 2
        xL1 = L(1,1); yL1 = L(1,2); zL1 = L(1,3);
        xL2 = L(2,1); yL2 = L(2,2); zL2 = L(2,3);
    elseif isnumeric(L) && numel(L) >= 2
        xL1 = L(1); yL1 = 0; zL1 = 0;
        xL2 = L(2); yL2 = 0; zL2 = 0;
    else
        error('Unsupported output format from astro.cr3bp.lagrangePoints.');
    end
end