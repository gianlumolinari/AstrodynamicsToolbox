clc;
clear;
close all;

% BUILD_HALO_FAMILY_DATABASE
% Build a small validated halo database from a continued family.

mu = 0.012150585609624;

%% ------------------------------------------------------------------------
% 1) Build family using two corrected JPL seeds
% -------------------------------------------------------------------------
libr = 1;
branchJpl = 'N';

jplIndex1 = 1;
jplIndex2 = 2;

nMembers = 20;
ds = 1e-3;

data = astro.cr3bp.queryJPLPeriodicOrbits( ...
    'sys', 'earth-moon', ...
    'family', 'halo', ...
    'libr', libr, ...
    'branch', branchJpl);

fields = string(data.fields);
rows = data.data;

seed1raw = localExtractSeed(rows{jplIndex1}, fields);
seed2raw = localExtractSeed(rows{jplIndex2}, fields);

corr1 = astro.cr3bp.differentialCorrectionHalo( ...
    seed1raw.x0, seed1raw.z0, seed1raw.vy0, 0.5*seed1raw.T, mu, 30, 1e-12);
corr2 = astro.cr3bp.differentialCorrectionHalo( ...
    seed2raw.x0, seed2raw.z0, seed2raw.vy0, 0.5*seed2raw.T, mu, 30, 1e-12);

if ~corr1.converged || ~corr2.converged
    error('Could not generate the initial corrected halo seeds.');
end

seed1 = struct('state', corr1.state0, 'period', corr1.period);
seed2 = struct('state', corr2.state0, 'period', corr2.period);

family = astro.cr3bp.continueHaloPseudoArc(seed1, seed2, nMembers, ds, mu);
nFam = numel(family);

fprintf('\nHalo family generated with %d members.\n', nFam);

%% ------------------------------------------------------------------------
% 2) Choose representative members
% -------------------------------------------------------------------------
idxSmall  = 1;
idxMedium = round((nFam+1)/2);
idxLarge  = nFam;

idxList = [idxSmall, idxMedium, idxLarge];
qualifiers = {'small', 'medium', 'large'};

%% ------------------------------------------------------------------------
% 3) Validate and save representative members
% -------------------------------------------------------------------------
for i = 1:numel(idxList)
    k = idxList(i);

    state0 = family(k).state0;
    T = family(k).period;

    orbit = astro.periodic.buildOrbitStruct(state0, T, mu, ...
        struct('family','halo', ...
               'libr', libr, ...
               'system','earth-moon', ...
               'source','halo family continuation representative member', ...
               'dimension','3D'), ...
        struct('verbose', false));

    report = astro.manifolds.validatePeriodicOrbitForManifolds(orbit, mu, struct('verbose', true));

    if ~report.isAcceptable
        warning('Family member %d failed validation. Skipping save.', k);
        continue
    end

    astro.periodic.saveValidatedOrbit(state0, T, mu, struct( ...
        'system', 'earth_moon', ...
        'family', 'halo', ...
        'libr', libr, ...
        'branch', 'north', ...
        'qualifier', qualifiers{i}, ...
        'source', sprintf('Halo family continuation member %d', k), ...
        'closureError', orbit.closureError));
end

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