clc;
clear;
close all;

%% CR3BP Periodic Orbits and Manifolds
% Test case: Sun-Earth planar Lyapunov family around L2

%% ------------------------------------------------------------------------
% User settings
% -------------------------------------------------------------------------
systemName = 'sun_earth';          % 'earth_moon', 'sun_earth'
libr       = 2;                    % 1 or 2
familyType = 'halo';    % 'planar_lyapunov', 'halo'
branch     = 'north';                   % planar Lyapunov: leave empty

seedSource = 'jpl';                % 'jpl' or 'validated'
validatedSeedName = '';            % only used if seedSource = 'validated'

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

%% ------------------------------------------------------------------------
% System constants
% -------------------------------------------------------------------------
mu = astro.cr3bp.getSystemMu(systemName);
L = astro.cr3bp.lagrangePoints(mu);

fprintf('\nConfiguration\n');
fprintf('-------------\n');
fprintf('System      : %s\n', systemName);
fprintf('Family      : %s\n', familyType);
fprintf('Libration   : L%d\n', libr);
fprintf('Branch      : %s\n', string(branch));
fprintf('Seed source : %s\n', seedSource);

%% ------------------------------------------------------------------------
% Build corrected seeds
% -------------------------------------------------------------------------
[seed1, seed2] = astro.periodic.buildSeedFromSource( ...
    systemName, familyType, libr, branch, seedSource, ...
    struct('validatedSeedName', validatedSeedName, ...
           'jplIndex1', jplIndex1, ...
           'jplIndex2', jplIndex2));

fprintf('\nSeed summary\n');
fprintf('------------\n');
fprintf('Seed 1 period = %.16f\n', seed1.period);
fprintf('Seed 2 period = %.16f\n', seed2.period);
fprintf('Seed 1 x0     = %.16f\n', seed1.state(1));
fprintf('Seed 2 x0     = %.16f\n', seed2.state(1));

% Sanity check: for L1/L2 planar Lyapunov, x0 should be near the chosen libration point
if libr == 1
    xL = L.L1(1);
elseif libr == 2
    xL = L.L2(1);
else
    error('Only L1/L2 supported here.');
end

fprintf('Target libration x = %.16f\n', xL);
fprintf('|seed1.x0 - xL|    = %.3e\n', abs(seed1.state(1) - xL));
fprintf('|seed2.x0 - xL|    = %.3e\n', abs(seed2.state(1) - xL));

%% ------------------------------------------------------------------------
% Validate first corrected orbit
% -------------------------------------------------------------------------
if doValidateSeedOrbit
    orbit0 = astro.periodic.buildOrbitStruct(seed1.state, seed1.period, mu, ...
        struct('family', familyType, ...
               'libr', libr, ...
               'system', strrep(systemName,'_','-'), ...
               'source', 'live workflow seed orbit', ...
               'dimension', localDimensionFromFamily(familyType)), ...
        struct('verbose', true));

    report0 = astro.manifolds.validatePeriodicOrbitForManifolds( ...
        orbit0, mu, struct('verbose', true));

    fprintf('\nSeed orbit validation summary\n');
    fprintf('-----------------------------\n');
    fprintf('Closure error : %.3e\n', orbit0.closureError);
    fprintf('Acceptable    : %d\n', report0.isAcceptable);
    fprintf('Seed orbit x-range: [%.6f, %.6f]\n', min(orbit0.x(:,1)), max(orbit0.x(:,1)));

    figure;
    plot(orbit0.x(:,1), orbit0.x(:,2), 'k', 'LineWidth', 2); hold on; grid on; axis equal;
    plot(-mu, 0, 'ko', 'MarkerFaceColor','k');
    plot(1-mu, 0, 'o', 'Color',[0.4 0.4 0.4], 'MarkerFaceColor',[0.7 0.7 0.7]);
    plot(L.L1(1), L.L1(2), 'yo', 'MarkerFaceColor','y');
    plot(L.L2(1), L.L2(2), 'co', 'MarkerFaceColor','c');
    xlabel('x'); ylabel('y');
    title(sprintf('Seed orbit check: %s, L%d', familyType, libr));

    % Zoom near selected libration point
    xlim([xL - 0.05, xL + 0.05]);
end

%% ------------------------------------------------------------------------
% Continue family
% -------------------------------------------------------------------------
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

%% ------------------------------------------------------------------------
% Save representative family members
% -------------------------------------------------------------------------
if doContinueFamily && doSaveRepresentativeMembers
    nFam = numel(family);

    idxSmall  = 1;
    idxMedium = round((nFam + 1)/2);
    idxLarge  = nFam;

    idxList = [idxSmall, idxMedium, idxLarge];
    qualifiers = {'small', 'medium', 'large'};

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
            warning('Representative member %d failed validation. Skipping save.', k);
            continue
        end

        fprintf('Saving %s: x-range = [%.6f, %.6f]\n', qualifiers{i}, ...
            min(orbitk.x(:,1)), max(orbitk.x(:,1)));

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

%% ------------------------------------------------------------------------
% Generate manifolds from a saved representative member
% -------------------------------------------------------------------------
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

    fprintf('Loaded representative x-range: [%.6f, %.6f]\n', ...
        min(orbit.x(:,1)), max(orbit.x(:,1)));

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
            fig = figure;
            ax = axes(fig);
            astro.plot.plotPeriodicOrbitManifolds2D( ...
                man, orbit, S.mu, sprintf('Live workflow manifolds: %s', orbitKey), ax);
            xlim(ax, [xL - 0.05, xL + 0.05]);

        case 'halo'
            astro.plot.plotPeriodicOrbitManifolds3D( ...
                man, orbit, S.mu, sprintf('Live workflow manifolds: %s', orbitKey));

        otherwise
            error('Unsupported family type for plotting: %s', familyType);
    end
end

%% ------------------------------------------------------------------------
% Local helper functions
% -------------------------------------------------------------------------
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