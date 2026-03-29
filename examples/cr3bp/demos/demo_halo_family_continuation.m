clc;
clear;
close all;

% DEMO_HALO_FAMILY_CONTINUATION
% Continue a halo family starting from two corrected JPL halo seeds.

mu = 0.012150585609624;

%% ------------------------------------------------------------------------
% 1) User settings
% -------------------------------------------------------------------------
libr = 1;
branch = 'N';

jplIndex1 = 1;
jplIndex2 = 2;

nMembers = 50;
ds = 1e-3;

%% ------------------------------------------------------------------------
% 2) Query JPL halo family
% -------------------------------------------------------------------------
data = astro.cr3bp.queryJPLPeriodicOrbits( ...
    'sys', 'earth-moon', ...
    'family', 'halo', ...
    'libr', libr, ...
    'branch', branch);

fields = string(data.fields);
rows = data.data;

if jplIndex1 > numel(rows) || jplIndex2 > numel(rows)
    error('Requested JPL indices exceed available halo rows.');
end

%% ------------------------------------------------------------------------
% 3) Build corrected seed 1
% -------------------------------------------------------------------------
seed1raw = localExtractSeed(rows{jplIndex1}, fields);

fprintf('JPL seed 1:\n');
fprintf('  row  = %d\n', jplIndex1);
fprintf('  x0   = %.16f\n', seed1raw.x0);
fprintf('  z0   = %.16f\n', seed1raw.z0);
fprintf('  vy0  = %.16f\n', seed1raw.vy0);
fprintf('  T    = %.16f\n', seed1raw.T);

corr1 = astro.cr3bp.differentialCorrectionHalo( ...
    seed1raw.x0, seed1raw.z0, seed1raw.vy0, 0.5*seed1raw.T, mu, 30, 1e-12);

if ~corr1.converged
    error('Failed to correct halo seed 1.');
end

seed1 = struct();
seed1.state = corr1.state0;
seed1.period = corr1.period;

fprintf('\nCorrected seed 1:\n');
fprintf('  x0   = %.16f\n', seed1.state(1));
fprintf('  z0   = %.16f\n', seed1.state(3));
fprintf('  vy0  = %.16f\n', seed1.state(5));
fprintf('  T    = %.16f\n', seed1.period);

%% ------------------------------------------------------------------------
% 4) Build corrected seed 2
% -------------------------------------------------------------------------
seed2raw = localExtractSeed(rows{jplIndex2}, fields);

fprintf('\nJPL seed 2:\n');
fprintf('  row  = %d\n', jplIndex2);
fprintf('  x0   = %.16f\n', seed2raw.x0);
fprintf('  z0   = %.16f\n', seed2raw.z0);
fprintf('  vy0  = %.16f\n', seed2raw.vy0);
fprintf('  T    = %.16f\n', seed2raw.T);

corr2 = astro.cr3bp.differentialCorrectionHalo( ...
    seed2raw.x0, seed2raw.z0, seed2raw.vy0, 0.5*seed2raw.T, mu, 30, 1e-12);

if ~corr2.converged
    error('Failed to correct halo seed 2.');
end

seed2 = struct();
seed2.state = corr2.state0;
seed2.period = corr2.period;

fprintf('\nCorrected seed 2:\n');
fprintf('  x0   = %.16f\n', seed2.state(1));
fprintf('  z0   = %.16f\n', seed2.state(3));
fprintf('  vy0  = %.16f\n', seed2.state(5));
fprintf('  T    = %.16f\n', seed2.period);

%% ------------------------------------------------------------------------
% 5) Continue family
% -------------------------------------------------------------------------
family = astro.cr3bp.continueHaloPseudoArc(seed1, seed2, nMembers, ds, mu);
nFam = numel(family);

fprintf('\nFamily continuation complete.\n');
fprintf('Stored members: %d\n', nFam);

%% ------------------------------------------------------------------------
% 6) Build orbit structs and diagnostics
% -------------------------------------------------------------------------
Tvals = zeros(nFam,1);
xvals = zeros(nFam,1);
zvals = zeros(nFam,1);
vyvals = zeros(nFam,1);
cerr = zeros(nFam,1);
Cvals = zeros(nFam,1);

for k = 1:nFam
    state0 = family(k).state0;
    T = family(k).period;

    orbitk = astro.periodic.buildOrbitStruct(state0, T, mu, ...
        struct('family','halo', ...
               'libr', libr, ...
               'system','earth-moon', ...
               'source','halo pseudo-arclength continuation', ...
               'dimension','3D'), ...
        struct('verbose', false));

    if k == 1
        orbits = repmat(orbitk, nFam, 1);
    end
    orbits(k) = orbitk;

    Tvals(k)  = T;
    xvals(k)  = state0(1);
    zvals(k)  = state0(3);
    vyvals(k) = state0(5);
    cerr(k)   = orbitk.closureError;
    Cvals(k)  = family(k).C;
end

fprintf('\nFamily diagnostics\n');
fprintf('------------------\n');
fprintf('Min closure error: %.3e\n', min(cerr));
fprintf('Max closure error: %.3e\n', max(cerr));
fprintf('Min period       : %.16f\n', min(Tvals));
fprintf('Max period       : %.16f\n', max(Tvals));
fprintf('Min Jacobi C     : %.16f\n', min(Cvals));
fprintf('Max Jacobi C     : %.16f\n', max(Cvals));

%% ------------------------------------------------------------------------
% 7) Plot family in 3D
% -------------------------------------------------------------------------
fig = figure;
ax = axes(fig);
hold(ax,'on');
grid(ax,'on');
box(ax,'on');

for k = 1:nFam
    X = orbits(k).x;
    plot3(ax, X(:,1), X(:,2), X(:,3), 'LineWidth', 2.0);
end

L = astro.cr3bp.lagrangePoints(mu);
plot3(ax, -mu, 0, 0, 'ko', 'MarkerFaceColor','k', 'MarkerSize', 6);
plot3(ax, 1-mu, 0, 0, 'o', ...
    'Color',[0.4 0.4 0.4], ...
    'MarkerFaceColor',[0.7 0.7 0.7], ...
    'MarkerSize', 8);
plot3(ax, L.L1(1), L.L1(2), L.L1(3), 'yo', 'MarkerFaceColor','y', 'MarkerSize', 7);
plot3(ax, L.L2(1), L.L2(2), L.L2(3), 'co', 'MarkerFaceColor','c', 'MarkerSize', 7);

axis(ax,'equal');
xlabel(ax,'x');
ylabel(ax,'y');
zlabel(ax,'z');
title(ax, sprintf('Halo family continuation from JPL seeds, L%d %s', libr, upper(branch)));
view(ax,[35 25]);

%% ------------------------------------------------------------------------
% 8) Trend plots
% -------------------------------------------------------------------------
figure;
plot(1:nFam, Tvals, 'o-', 'LineWidth', 1.5);
grid on;
xlabel('Family member index');
ylabel('Period');
title('Halo family period trend');

figure;
plot(1:nFam, zvals, 'o-', 'LineWidth', 1.5);
grid on;
xlabel('Family member index');
ylabel('z_0');
title('Halo family z_0 trend');

figure;
plot(1:nFam, xvals, 'o-', 'LineWidth', 1.5);
grid on;
xlabel('Family member index');
ylabel('x_0');
title('Halo family x_0 trend');

figure;
plot(1:nFam, Cvals, 'o-', 'LineWidth', 1.5);
grid on;
xlabel('Family member index');
ylabel('Jacobi constant');
title('Halo family Jacobi constant trend');

figure;
semilogy(1:nFam, cerr, 'o-', 'LineWidth', 1.5);
grid on;
xlabel('Family member index');
ylabel('Closure error');
title('Halo family closure error');

%% local helpers
function seed = localExtractSeed(rawRow, fields)
    row = rawRow;
    while iscell(row) && numel(row) == 1
        row = row{1};
    end
    row = row(:).';

    seed = struct();
    seed.x0  = localGetFieldValue(row, fields, 'x');
    seed.y0  = localGetFieldValue(row, fields, 'y');
    seed.z0  = localGetFieldValue(row, fields, 'z');
    seed.vx0 = localGetFieldValue(row, fields, 'vx');
    seed.vy0 = localGetFieldValue(row, fields, 'vy');
    seed.vz0 = localGetFieldValue(row, fields, 'vz');
    seed.T   = localGetFieldValue(row, fields, 'period');
    seed.C   = localGetFieldValue(row, fields, 'jacobi');
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