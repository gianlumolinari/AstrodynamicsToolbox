clc;
clear;
close all;

% BUILD_PLANAR_LYAPUNOV_FAMILY_DATABASE
% Build a small validated planar Lyapunov database from a continued family.

mu = 0.012150585609624;
libr = 1;

%% ------------------------------------------------------------------------
% 1) Load a validated planar Lyapunov orbit as starting seed
% -------------------------------------------------------------------------
S = astro.periodic.loadValidatedOrbit('planar_lyapunov_L1');

x0 = S.state0(1);
vy0 = S.state0(5);
thalf = 0.5 * S.T;

seed1 = struct();
seed1.state = S.state0;
seed1.period = S.T;

fprintf('Validated planar seed:\n');
fprintf('  x0    = %.16f\n', x0);
fprintf('  vy0   = %.16f\n', vy0);
fprintf('  T/2   = %.16f\n', thalf);
fprintf('  T     = %.16f\n', S.T);

%% ------------------------------------------------------------------------
% 2) Generate second nearby corrected seed
% -------------------------------------------------------------------------
uPred2 = [x0 + 1e-3; vy0; thalf];
tangent0 = [1; 0; 0];
tangent0 = tangent0 / norm(tangent0);

corr2 = astro.cr3bp.correctPlanarLyapunovPseudoArc(uPred2, tangent0, mu, 20, 1e-11);

if ~corr2.converged
    error('Failed to generate second corrected planar Lyapunov seed.');
end

seed2 = struct();
seed2.state = corr2.state0;
seed2.period = corr2.period;

fprintf('\nSecond corrected planar seed:\n');
fprintf('  x0    = %.16f\n', corr2.u(1));
fprintf('  vy0   = %.16f\n', corr2.u(2));
fprintf('  T/2   = %.16f\n', corr2.u(3));
fprintf('  T     = %.16f\n', corr2.period);

%% ------------------------------------------------------------------------
% 3) Continue the planar family
% -------------------------------------------------------------------------
nMembers = 20;
ds = 1e-3;

family = astro.cr3bp.continueFamilyPseudoArc(seed1, seed2, nMembers, ds, mu);
nFam = numel(family);

fprintf('\nPlanar family generated with %d members.\n', nFam);

%% ------------------------------------------------------------------------
% 4) Choose representative members
% -------------------------------------------------------------------------
idxSmall  = 1;
idxMedium = round((nFam+1)/2);
idxLarge  = nFam;

idxList = [idxSmall, idxMedium, idxLarge];
qualifiers = {'small', 'medium', 'large'};

%% ------------------------------------------------------------------------
% 5) Validate and save representative members
% -------------------------------------------------------------------------
for i = 1:numel(idxList)
    k = idxList(i);

    state0 = family(k).state0;
    T = family(k).period;

    orbit = astro.periodic.buildOrbitStruct(state0, T, mu, ...
        struct('family','planar_lyapunov', ...
               'libr', libr, ...
               'system','earth-moon', ...
               'source','planar Lyapunov family continuation representative member', ...
               'dimension','2D'), ...
        struct('verbose', false));

    report = astro.manifolds.validatePeriodicOrbitForManifolds(orbit, mu, struct('verbose', true));

    if ~report.isAcceptable
        warning('Planar family member %d failed validation. Skipping save.', k);
        continue
    end

    astro.periodic.saveValidatedOrbit(state0, T, mu, struct( ...
        'system', 'earth_moon', ...
        'family', 'planar_lyapunov', ...
        'libr', libr, ...
        'branch', '', ...
        'qualifier', qualifiers{i}, ...
        'source', sprintf('Planar Lyapunov family continuation member %d', k), ...
        'closureError', orbit.closureError));
end