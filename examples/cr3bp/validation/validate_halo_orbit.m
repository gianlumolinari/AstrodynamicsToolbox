clc;
clear;
close all;

% VALIDATE_HALO_ORBIT
% Validate one halo orbit from a JPL seed + local differential correction,
% then save it as a validated orbit if closure is good.

mu = 0.012150585609624;

%% USER OPTIONS
libr = 1;            % 1 or 2
branch = 'N';        % 'N' or 'S'
jplIndex = 100;        % candidate row from JPL query
saveOrbit = true;

%% ------------------------------------------------------------------------
% 1) Query JPL halo seed
% -------------------------------------------------------------------------
data = astro.cr3bp.queryJPLPeriodicOrbits( ...
    'sys', 'earth-moon', ...
    'family', 'halo', ...
    'libr', libr, ...
    'branch', branch);

fields = string(data.fields);
rows = data.data;

if jplIndex > numel(rows)
    error('Requested JPL row %d but only %d rows are available.', jplIndex, numel(rows));
end

row = localNormalizeRow(rows{jplIndex});

seed = struct();
seed.x0  = localGetFieldValue(row, fields, 'x');
seed.y0  = localGetFieldValue(row, fields, 'y');
seed.z0  = localGetFieldValue(row, fields, 'z');
seed.vx0 = localGetFieldValue(row, fields, 'vx');
seed.vy0 = localGetFieldValue(row, fields, 'vy');
seed.vz0 = localGetFieldValue(row, fields, 'vz');
seed.T   = localGetFieldValue(row, fields, 'period');
seed.C   = localGetFieldValue(row, fields, 'jacobi');

fprintf('Using JPL halo seed:\n');
fprintf('  libr   = L%d\n', libr);
fprintf('  branch = %s\n', branch);
fprintf('  row    = %d\n', jplIndex);
fprintf('  x0     = %.16f\n', seed.x0);
fprintf('  y0     = %.16f\n', seed.y0);
fprintf('  z0     = %.16f\n', seed.z0);
fprintf('  vx0    = %.16f\n', seed.vx0);
fprintf('  vy0    = %.16f\n', seed.vy0);
fprintf('  vz0    = %.16f\n', seed.vz0);
fprintf('  T      = %.16f\n', seed.T);
fprintf('  C      = %.16f\n', seed.C);

%% ------------------------------------------------------------------------
% 2) Build symmetry-compatible initial guess
% -------------------------------------------------------------------------
% Your corrector assumes:
%   X0 = [x0; 0; z0; 0; vy0; 0]
% and corrects [x0, vy0, thalf], keeping z0 fixed.
%
% So we use the JPL row as a seed, but project it to the symmetry form.

x0Guess    = seed.x0;
z0Guess    = seed.z0;
vy0Guess   = seed.vy0;
thalfGuess = 0.5 * seed.T;

fprintf('\nProjected symmetry-form guess:\n');
fprintf('  x0_guess    = %.16f\n', x0Guess);
fprintf('  z0_guess    = %.16f\n', z0Guess);
fprintf('  vy0_guess   = %.16f\n', vy0Guess);
fprintf('  thalf_guess = %.16f\n', thalfGuess);

%% ------------------------------------------------------------------------
% 3) Halo differential correction
% -------------------------------------------------------------------------
corr = astro.cr3bp.differentialCorrectionHalo( ...
    x0Guess, z0Guess, vy0Guess, thalfGuess, mu, 30, 1e-12);

if ~corr.converged
    error('Halo differential correction did not converge.');
end

state0 = corr.state0;
T = corr.period;

fprintf('\nCorrected halo orbit\n');
fprintf('--------------------\n');
fprintf('x0  = %.16f\n', state0(1));
fprintf('y0  = %.16f\n', state0(2));
fprintf('z0  = %.16f\n', state0(3));
fprintf('vx0 = %.16f\n', state0(4));
fprintf('vy0 = %.16f\n', state0(5));
fprintf('vz0 = %.16f\n', state0(6));
fprintf('T   = %.16f\n', T);

orbit = astro.periodic.buildOrbitStruct(state0, T, mu, ...
    struct('family','halo', ...
           'libr', libr, ...
           'system','earth-moon', ...
           'source', sprintf('JPL halo seed row %d corrected locally', jplIndex), ...
           'dimension','3D'), ...
    struct('verbose', true));

%% ------------------------------------------------------------------------
% 4) Validation
% -------------------------------------------------------------------------
report = astro.manifolds.validatePeriodicOrbitForManifolds(orbit, mu, struct('verbose', true));

xf = orbit.x(end,:).';
dx = xf - orbit.state0;

fprintf('State closure mismatch:\n');
fprintf('  dx  = %.3e\n', dx(1));
fprintf('  dy  = %.3e\n', dx(2));
fprintf('  dz  = %.3e\n', dx(3));
fprintf('  dvx = %.3e\n', dx(4));
fprintf('  dvy = %.3e\n', dx(5));
fprintf('  dvz = %.3e\n', dx(6));

%% ------------------------------------------------------------------------
% 5) Save if validated
% -------------------------------------------------------------------------
if saveOrbit && report.isAcceptable
    if libr == 1 && strcmpi(branch,'N')
        outFile = 'data/validated/cr3bp/earth_moon/halo_L1_north_validated.mat';
    elseif libr == 1 && strcmpi(branch,'S')
        outFile = 'data/validated/cr3bp/earth_moon/halo_L1_south_validated.mat';
    elseif libr == 2 && strcmpi(branch,'N')
        outFile = 'data/validated/cr3bp/earth_moon/halo_L2_north_validated.mat';
    elseif libr == 2 && strcmpi(branch,'S')
        outFile = 'data/validated/cr3bp/earth_moon/halo_L2_south_validated.mat';
    else
        outFile = sprintf('data/validated/cr3bp/earth_moon/halo_L%d_%s_validated.mat', ...
            libr, lower(branch));
    end

    family = 'halo';
    closureError = orbit.closureError;
    source = orbit.source;

    save(outFile, 'state0', 'T', 'mu', 'family', 'libr', 'branch', 'source', 'closureError');
    fprintf('Saved validated halo orbit to %s\n', outFile);
else
    fprintf('Halo orbit not saved.\n');
end

%% ------------------------------------------------------------------------
% 6) Plot corrected orbit
% -------------------------------------------------------------------------
fig = figure;
ax = axes(fig);
hold(ax,'on');
grid(ax,'on');
box(ax,'on');

plot3(ax, orbit.x(:,1), orbit.x(:,2), orbit.x(:,3), 'k', 'LineWidth', 2.0);

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
title(ax, sprintf('Validated halo orbit: L%d %s', libr, upper(branch)));
view(ax, [35 25]);

%% local helpers
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