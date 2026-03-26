clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/astroToolbox')
startup

% ==========================================================
% JUICE full heliocentric comparison:
% piecewise HF propagation between gravity-assist anchors
% vs full SPICE truth
% ==========================================================

% Reset SPICE and load kernels
cspice_kclear;
astro.ephem.loadSpiceKernels(fullfile(pwd, 'data', 'spice'));

juiceKernelDir = fullfile(pwd, 'data', 'spice', 'juice');
if exist(juiceKernelDir, 'dir') ~= 7
    error('JUICE mission kernel directory not found:\n%s', juiceKernelDir);
end
astro.ephem.loadSpiceKernels(juiceKernelDir);

% ----------------------------------------------------------
% Bodies
% ----------------------------------------------------------
mercury = astro.bodies.getBody('mercury');
venus   = astro.bodies.getBody('venus');
earth   = astro.bodies.getBody('earth');
moon    = astro.bodies.getBody('moon');
mars    = astro.bodies.getBody('mars');
jupiter = astro.bodies.getBody('jupiter');
saturn  = astro.bodies.getBody('saturn');

% ----------------------------------------------------------
% Resolve JUICE target
% ----------------------------------------------------------
candidateTargets = {'JUICE', 'JUICE SPACECRAFT', '-28'};
juiceTarget = '';

for k = 1:numel(candidateTargets)
    try
        tmp = astro.ephem.getSpiceState(candidateTargets{k}, ...
            '2023-04-14 12:14:00', 'SUN', 'J2000', 'NONE');
        juiceTarget = tmp.target;
        break
    catch
    end
end

if isempty(juiceTarget)
    error('Could not resolve JUICE spacecraft target name automatically.');
end

fprintf('Resolved JUICE spacecraft target as: %s\n', juiceTarget);

% ----------------------------------------------------------
% Segment anchors
% These are the epochs at which the trajectory is reset from SPICE,
% so the gravity assists are incorporated segment by segment.
% ----------------------------------------------------------
epochs = { ...
    '2023-04-14 12:14:00'   % launch
    '2024-08-20 21:56:00'   % LEGA anchor
    '2025-08-31 00:00:00'   % Venus flyby anchor
    '2026-09-01 00:00:00'   % Earth flyby anchor
    '2029-01-01 00:00:00'   % Earth flyby anchor
    '2031-07-01 00:00:00'   % Jupiter arrival anchor
    };

segmentNames = { ...
    'Launch to LEGA'
    'LEGA to Venus'
    'Venus to Earth'
    'Earth to Earth'
    'Earth to Jupiter'
    };

% ----------------------------------------------------------
% High-fidelity propagation config
% ----------------------------------------------------------
configHF = struct();
configHF.RelTol   = 1e-12;
configHF.AbsTol   = 1e-12;
configHF.Cr       = 1.2;
configHF.A_over_m = 0.0;

configHF.bodyNames = { ...
    'MERCURY BARYCENTER', ...
    'VENUS BARYCENTER', ...
    'EARTH BARYCENTER', ...
    'MOON', ...
    'MARS BARYCENTER', ...
    'JUPITER BARYCENTER', ...
    'SATURN BARYCENTER'};

configHF.muBodies = [ ...
    mercury.mu, ...
    venus.mu, ...
    earth.mu, ...
    moon.mu, ...
    mars.mu, ...
    jupiter.mu, ...
    saturn.mu];

% ----------------------------------------------------------
% Piecewise HF propagation
% ----------------------------------------------------------
fprintf('\n============================================================\n');
fprintf('Piecewise high-fidelity propagation through assist anchors\n');
fprintf('============================================================\n');

xHFfull = [];
tHFfull = [];

for i = 1:(numel(epochs)-1)

    t0UTC = epochs{i};
    tfUTC = epochs{i+1};

    et0 = cspice_str2et(t0UTC);
    etf = cspice_str2et(tfUTC);
    tof = etf - et0;

    sc0 = astro.ephem.getSpiceState(juiceTarget, t0UTC, 'SUN', 'J2000', 'NONE');
    x0 = [sc0.r; sc0.v];

    fprintf('\nSegment %d: %s\n', i, segmentNames{i});
    fprintf('  Start: %s\n', t0UTC);
    fprintf('  End  : %s\n', tfUTC);
    fprintf('  TOF  : %.2f days\n', tof/86400);

    outHF = astro.propagators.propagateHighFidelity(x0, et0, tof, configHF);

    if i == 1
        xHFfull = outHF.x;
        tHFfull = et0 + outHF.t;
    else
        xHFfull = [xHFfull; outHF.x(2:end,:)]; %#ok<AGROW>
        tHFfull = [tHFfull; et0 + outHF.t(2:end)]; %#ok<AGROW>
    end

    scf = astro.ephem.getSpiceState(juiceTarget, tfUTC, 'SUN', 'J2000', 'NONE');
    xEnd = outHF.x(end,:).';

    posErr = norm(xEnd(1:3) - scf.r);
    velErr = norm(xEnd(4:6) - scf.v);

    fprintf('  Terminal position error : %12.3f km\n', posErr);
    fprintf('  Terminal velocity error : %12.6f km/s\n', velErr);
end

% ----------------------------------------------------------
% Full SPICE truth sampled over same stitched epochs
% ----------------------------------------------------------
fprintf('\nSampling SPICE truth along stitched propagated epochs...\n');

xTruthFull = zeros(numel(tHFfull), 6);

for k = 1:numel(tHFfull)
    utcK = cspice_et2utc(tHFfull(k), 'C', 0);
    scK = astro.ephem.getSpiceState(juiceTarget, utcK, 'SUN', 'J2000', 'NONE');
    xTruthFull(k,:) = [scK.r(:).' scK.v(:).'];
end

% ----------------------------------------------------------
% Global errors along stitched trajectory
% ----------------------------------------------------------
posErrHist = vecnorm(xHFfull(:,1:3) - xTruthFull(:,1:3), 2, 2);
velErrHist = vecnorm(xHFfull(:,4:6) - xTruthFull(:,4:6), 2, 2);

fprintf('\n============================================================\n');
fprintf('Global stitched-trajectory comparison\n');
fprintf('============================================================\n');
fprintf('  Final position error : %12.3f km\n', posErrHist(end));
fprintf('  Final velocity error : %12.6f km/s\n', velErrHist(end));
fprintf('  Max position error   : %12.3f km\n', max(posErrHist));
fprintf('  Max velocity error   : %12.6f km/s\n', max(velErrHist));
fprintf('  RMS position error   : %12.3f km\n', sqrt(mean(posErrHist.^2)));
fprintf('  RMS velocity error   : %12.6f km/s\n', sqrt(mean(velErrHist.^2)));

% ----------------------------------------------------------
% Plot full heliocentric comparison
% ----------------------------------------------------------
figure('Color','w');
hold on;
grid on;
axis equal;

plot(0, 0, 'ko', ...
    'MarkerSize', 10, ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Sun');

plot(xTruthFull(:,1), xTruthFull(:,2), 'k--', ...
    'LineWidth', 1.6, ...
    'DisplayName', 'JUICE SPICE truth');

plot(xHFfull(:,1), xHFfull(:,2), 'b-', ...
    'LineWidth', 1.8, ...
    'DisplayName', 'Piecewise HF propagated');

% Anchor markers
for i = 1:numel(epochs)
    scA = astro.ephem.getSpiceState(juiceTarget, epochs{i}, 'SUN', 'J2000', 'NONE');
    plot(scA.r(1), scA.r(2), 'ro', ...
        'MarkerSize', 6, ...
        'LineWidth', 1.2, ...
        'HandleVisibility', 'off');
end

xlabel('x [km]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('y [km]', 'FontSize', 13, 'FontWeight', 'bold');
title('JUICE Full Heliocentric Trajectory: Piecewise HF vs SPICE', ...
    'FontSize', 15, 'FontWeight', 'bold');

legend('Location', 'best');

fprintf('\nDone.\n');