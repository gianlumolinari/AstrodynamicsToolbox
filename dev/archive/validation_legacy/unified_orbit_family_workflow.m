clc;
clear;
close all;

% UNIFIED_ORBIT_FAMILY_WORKFLOW
% Generic workflow for:
%   1) selecting system / libration point / family
%   2) building corrected seeds
%   3) validating a seed orbit
%   4) continuing a family
%   5) saving representative members
%   6) generating manifolds from a saved representative member
%
% Supported families for now:
%   - planar_lyapunov
%   - halo
%
% Supported systems for now:
%   - earth_moon
%   - sun_earth   (subject to seed/catalog availability)

%% ========================================================================
% USER SETTINGS
% ========================================================================
systemName = 'earth_moon';          % 'earth_moon', 'sun_earth'
libr       = 1;                     % 1 or 2
familyType = 'halo';                % 'planar_lyapunov', 'halo'
branch     = 'north';               % '', 'north', 'south'

seedSource = 'jpl';                 % 'jpl' or 'validated'
validatedSeedName = 'halo_L1_north';

jplIndex1 = 1;
jplIndex2 = 2;

doValidateSeedOrbit         = true;
doContinueFamily            = true;
doSaveRepresentativeMembers = true;
doGenerateManifolds         = true;

nMembers = 20;
ds = 1e-3;

eps0 = 2e-5;
tfStable = 3.0;
tfUnstable = 3.0;
numPhaseSamples = 30;

representativeQualifier = 'medium'; % 'validated', 'small', 'medium', 'large'

%% ========================================================================
% 1) SYSTEM CONSTANTS
% ========================================================================
mu = astro.cr3bp.getSystemMu(systemName);

fprintf('\nUnified workflow configuration\n');
fprintf('------------------------------\n');
fprintf('System      : %s\n', systemName);
fprintf('Family      : %s\n', familyType);
fprintf('Libration   : L%d\n', libr);
fprintf('Branch      : %s\n', string(branch));
fprintf('Seed source : %s\n', seedSource);

%% ========================================================================
% 2) BUILD CORRECTED SEEDS
% ========================================================================
[seed1, seed2] = astro.periodic.buildSeedFromSource( ...
    systemName, familyType, libr, branch, seedSource, ...
    struct('validatedSeedName', validatedSeedName, ...
           'jplIndex1', jplIndex1, ...
           'jplIndex2', jplIndex2));

fprintf('\nSeed summary\n');
fprintf('------------\n');
fprintf('Seed 1 period = %.16f\n', seed1.period);
fprintf('Seed 2 period = %.16f\n', seed2.period);

%% ========================================================================
% 3) VALIDATE FIRST CORRECTED ORBIT
% ========================================================================
if doValidateSeedOrbit
    orbit0 = astro.periodic.buildOrbitStruct(seed1.state, seed1.period, mu, ...
        struct('family', familyType, ...
               'libr', libr, ...
               'system', strrep(systemName,'_','-'), ...
               'source', 'unified workflow seed orbit', ...
               'dimension', localDimensionFromFamily(familyType)), ...
        struct('verbose', true));

    report0 = astro.manifolds.validatePeriodicOrbitForManifolds( ...
        orbit0, mu, struct('verbose', true));

    fprintf('\nSeed orbit validation summary\n');
    fprintf('-----------------------------\n');
    fprintf('Closure error : %.3e\n', orbit0.closureError);
    fprintf('Acceptable    : %d\n', report0.isAcceptable);
end

%% ========================================================================
% 4) CONTINUE FAMILY
% ========================================================================
family = struct([]);

if doContinueFamily
    switch lower(familyType)
        case 'planar_lyapunov'
            family = astro.cr3bp.continueFamilyPseudoArc(seed1, seed2, nMembers, ds, mu);

        case 'halo'
            family = astro.cr3bp.continueHaloPseudoArc(seed1, seed2, nMembers, ds, mu);

        otherwise
            error('Unsupported family type for continuation: %s', familyType);
    end

    fprintf('\nFamily continuation complete.\n');
    fprintf('Stored members: %d\n', numel(family));
end

%% ========================================================================
% 5) SAVE REPRESENTATIVE FAMILY MEMBERS
% ========================================================================
if doContinueFamily && doSaveRepresentativeMembers
    nFam = numel(family);

    idxSmall  = 1;
    idxMedium = round((nFam + 1)/2);
    idxLarge  = nFam;

    idxList = [idxSmall, idxMedium, idxLarge];
    qualifiers = {'small', 'medium', 'large'};

    fprintf('\nSaving representative family members\n');
    fprintf('------------------------------------\n');

    for i = 1:numel(idxList)
        k = idxList(i);

        state0 = family(k).state0;
        T = family(k).period;

        orbitk = astro.periodic.buildOrbitStruct(state0, T, mu, ...
            struct('family', familyType, ...
                   'libr', libr, ...
                   'system', strrep(systemName,'_','-'), ...
                   'source', sprintf('%s family continuation member %d', familyType, k), ...
                   'dimension', localDimensionFromFamily(familyType)), ...
            struct('verbose', false));

        reportk = astro.manifolds.validatePeriodicOrbitForManifolds( ...
            orbitk, mu, struct('verbose', true));

        if ~reportk.isAcceptable
            warning('Representative family member %d failed validation. Skipping save.', k);
            continue
        end

        astro.periodic.saveValidatedOrbit(state0, T, mu, struct( ...
            'system', systemName, ...
            'family', familyType, ...
            'libr', libr, ...
            'branch', localBranchForSave(familyType, branch), ...
            'qualifier', qualifiers{i}, ...
            'source', sprintf('%s family continuation member %d', familyType, k), ...
            'closureError', orbitk.closureError));
    end
end

%% ========================================================================
% 6) LOAD ONE SAVED REPRESENTATIVE AND GENERATE MANIFOLDS
% ========================================================================
if doGenerateManifolds
    orbitKey = astro.periodic.makeOrbitKey( ...
        familyType, libr, localBranchForSave(familyType, branch), representativeQualifier);

    fprintf('\nGenerating manifolds from saved orbit: %s\n', orbitKey);

    S = astro.periodic.loadValidatedOrbit(orbitKey, systemName);
    orbit = astro.periodic.buildOrbitStruct(S.state0, S.T, S.mu, ...
        struct('family', S.family, ...
               'libr', S.libr, ...
               'system', strrep(systemName,'_','-'), ...
               'source', S.source, ...
               'dimension', localDimensionFromFamily(familyType)), ...
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

    switch lower(familyType)
        case 'planar_lyapunov'
            astro.plot.plotPeriodicOrbitManifolds2D( ...
                man, orbit, S.mu, sprintf('Unified workflow manifolds: %s', orbitKey));

        case 'halo'
            astro.plot.plotPeriodicOrbitManifolds3D( ...
                man, orbit, S.mu, sprintf('Unified workflow manifolds: %s', orbitKey));

        otherwise
            error('Unsupported family type for plotting: %s', familyType);
    end
end

%% ========================================================================
% LOCAL HELPERS
% ========================================================================
function dim = localDimensionFromFamily(familyType)
    switch lower(familyType)
        case 'planar_lyapunov'
            dim = '2D';
        case 'halo'
            dim = '3D';
        otherwise
            dim = '3D';
    end
end

function b = localBranchForSave(familyType, branch)
    if strcmpi(familyType, 'halo')
        b = lower(string(branch));
    else
        b = '';
    end
end