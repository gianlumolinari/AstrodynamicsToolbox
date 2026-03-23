clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/AstroToolbox')
startup
astro.ephem.loadSpiceKernels(fullfile(pwd, 'data', 'spice'));

% ==========================================================
% JUICE reference trajectory reconstruction
%
% Heliocentric MGA-style reconstruction of the published JUICE
% gravity-assist chronology.
%
% Sequence used here:
%   Earth launch         : 2023-04-14
%   Earth flyby (LEGA)   : 2024-08-20
%   Venus flyby          : 2025-08-31
%   Earth flyby          : 2026-09-01
%   Earth flyby          : 2029-01-01
%   Jupiter arrival      : 2031-07-01
%
% IMPORTANT:
%   The August 2024 lunar-Earth gravity assist is modeled here as a
%   single Earth-system flyby event in the heliocentric chain.
%   We do NOT insert a separate Moon->Earth heliocentric Lambert leg,
%   because that is not appropriate for a Sun-centered patched-conics
%   reconstruction and leads to a degenerate Lambert problem.
% ==========================================================

sun     = astro.bodies.getBody('sun');
earth   = astro.bodies.getBody('earth');
venus   = astro.bodies.getBody('venus');
jupiter = astro.bodies.getBody('jupiter');

% ----------------------------------------------------------
% Heliocentric sequence
% ----------------------------------------------------------
sequence = {'earth','earth','venus','earth','earth','jupiter'};

dates = {
    '2023-04-14 12:14:00'   % launch
    '2024-08-20 21:56:00'   % LEGA Earth flyby anchor
    '2025-08-31 00:00:00'   % Venus flyby
    '2026-09-01 00:00:00'   % Earth flyby anchor
    '2029-01-01 00:00:00'   % Earth flyby anchor
    '2031-07-01 00:00:00'   % Jupiter arrival anchor
};

fprintf('\n');
fprintf('============================================================\n');
fprintf('JUICE Reference Trajectory Reconstruction\n');
fprintf('============================================================\n');

fprintf('\nEncounter epochs used:\n');
for k = 1:numel(sequence)
    fprintf('  %-8s : %s\n', upper(sequence{k}), dates{k});
end

fprintf('\nNote: the 2024 lunar-Earth gravity assist is represented here\n');
fprintf('as a single Earth-system flyby in the heliocentric chain.\n');

% ----------------------------------------------------------
% Propagate heliocentric Lambert legs
% ----------------------------------------------------------
traj = astro.mga.propagateMGA(sequence, dates, sun.mu, 'spice');

% ----------------------------------------------------------
% Flyby feasibility setup
% ----------------------------------------------------------
bodyMap.earth = earth;
bodyMap.venus = venus;

minAltMap.earth = 300;   % km
minAltMap.venus = 300;   % km

report = astro.mga.evaluateMGA(traj, bodyMap, minAltMap, 1e-3);

% ----------------------------------------------------------
% Leg summary
% ----------------------------------------------------------
fprintf('\nLeg summary:\n');
for k = 1:numel(traj.legs)
    fprintf('\nLeg %d: %s -> %s\n', ...
        k, upper(traj.legs(k).bodyDepart), upper(traj.legs(k).bodyArrive));
    fprintf('  TOF                 : %.2f days\n', traj.legs(k).tofDays);
    fprintf('  |v_inf,dep|         : %.6f km/s\n', traj.legs(k).vInfDepMag);
    fprintf('  |v_inf,arr|         : %.6f km/s\n', traj.legs(k).vInfArrMag);
end

fprintf('\nFlyby diagnostic summary:\n');
for k = 1:report.nFlybys
    fb = report.flybys(k);
    chk = fb.check;

    fprintf('\nFlyby at %s\n', upper(fb.body));
    fprintf('  Incoming |v_inf^-|       : %.6f km/s\n', chk.magIn);
    fprintf('  Required |v_inf^+|       : %.6f km/s\n', chk.magOutReq);
    fprintf('  Magnitude mismatch       : %.6f km/s\n', chk.magMismatch);
    fprintf('  Required turn angle      : %.6f deg\n', chk.deltaReqDeg);
    fprintf('  Max turn angle @ rpMin   : %.6f deg\n', chk.deltaMaxDeg);
    fprintf('  Feasible (ballistic)     : %d\n', chk.feasible);
end

% ----------------------------------------------------------
% Heliocentric overview plot
% ----------------------------------------------------------
figure('Color','w');
hold on
axis equal
grid on

plot(0, 0, 'o', 'MarkerSize', 10, 'LineWidth', 1.5);

for k = 1:numel(traj.legs)
    x0Leg = [traj.legs(k).r1; traj.legs(k).v1];
    tspan = [0, traj.legs(k).tofSec];

    optsProp.RelTol = 1e-11;
    optsProp.AbsTol = 1e-11;
    optsProp.Solver = 'ode113';

    out = astro.propagators.propagate( ...
        @(t,x) astro.propagators.eomTwoBody(t, x, sun.mu), ...
        tspan, x0Leg, optsProp);

    astro.plot.plotOrbit2D(out.x(:,1:3), 'LineWidth', 1.3);
    plot(traj.legs(k).r1(1), traj.legs(k).r1(2), 'o', 'MarkerSize', 7, 'LineWidth', 1.2);
    plot(traj.legs(k).r2(1), traj.legs(k).r2(2), 's', 'MarkerSize', 7, 'LineWidth', 1.2);
end

xlabel('x [km]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('y [km]', 'FontSize', 13, 'FontWeight', 'bold');
title('JUICE Reference Trajectory: Heliocentric Overview', ...
    'FontSize', 15, 'FontWeight', 'bold');
legend('Sun', 'Leg arcs / encounters', 'Location', 'best');

% ----------------------------------------------------------
% Flyby plots
% ----------------------------------------------------------
rpEarth = earth.radius + minAltMap.earth;
rpVenus = venus.radius + minAltMap.venus;

% LEGA Earth-system flyby
astro.plot.plotFlyby2D( ...
    traj.legs(1).vInfArr, ...
    traj.legs(2).vInfDep, ...
    earth.mu, rpEarth, +1, earth.radius, 'Earth (LEGA)');

% Venus flyby
astro.plot.plotFlyby2D( ...
    traj.legs(2).vInfArr, ...
    traj.legs(3).vInfDep, ...
    venus.mu, rpVenus, +1, venus.radius, 'Venus');

% Earth flyby
astro.plot.plotFlyby2D( ...
    traj.legs(3).vInfArr, ...
    traj.legs(4).vInfDep, ...
    earth.mu, rpEarth, +1, earth.radius, 'Earth');

% Earth flyby
astro.plot.plotFlyby2D( ...
    traj.legs(4).vInfArr, ...
    traj.legs(5).vInfDep, ...
    earth.mu, rpEarth, +1, earth.radius, 'Earth');

% ----------------------------------------------------------
% Mission-level summary
% ----------------------------------------------------------
launchC3 = traj.legs(1).vInfDepMag^2;
arrivalVInfJupiter = traj.legs(end).vInfArrMag;
totalTOFdays = sum([traj.legs.tofDays]);

fprintf('\nMission-level summary:\n');
fprintf('  Launch C3                 : %.6f km^2/s^2\n', launchC3);
fprintf('  Jupiter arrival v_inf     : %.6f km/s\n', arrivalVInfJupiter);
fprintf('  Total time of flight      : %.2f days\n', totalTOFdays);
fprintf('  Total time of flight      : %.2f years\n', totalTOFdays / 365.25);

% ----------------------------------------------------------
% Optional validation against JUICE spacecraft SPK
% ----------------------------------------------------------
juiceKernelDir = fullfile(pwd, 'data', 'spice', 'juice');
hasJuiceDir = exist(juiceKernelDir, 'dir') == 7;

if hasJuiceDir
    fprintf('\nJUICE validation block enabled.\n');
    fprintf('Attempting to load JUICE mission kernels from:\n  %s\n', juiceKernelDir);

    try
        localLoadJuiceKernels(juiceKernelDir);

        fprintf('\nReconstruction vs JUICE SPK at encounter epochs:\n');
        fprintf('------------------------------------------------\n');

        for k = 1:numel(dates)
            epochUTC = dates{k};

            scState = [];
            targetNames = {'JUICE', 'JUICE SPACECRAFT', '-28'};
            for iTry = 1:numel(targetNames)
                try
                    scState = astro.ephem.getSpiceState(targetNames{iTry}, epochUTC, 'SUN', 'J2000', 'NONE');
                    break
                catch
                end
            end

            if isempty(scState)
                fprintf('  %s : could not resolve JUICE target name in loaded kernels.\n', epochUTC);
                continue
            end

            if k == 1
                rRef = traj.legs(1).r1;
            elseif k == numel(dates)
                rRef = traj.legs(end).r2;
            else
                rRef = traj.legs(k-1).r2;
            end

            posErr = norm(scState.r - rRef);

            fprintf('  %-19s  position difference = %.3f km\n', epochUTC, posErr);
        end

    catch ME
        fprintf('\nJUICE validation block could not be completed.\n');
        fprintf('Reason: %s\n', ME.message);
    end
else
    fprintf('\nNo JUICE mission SPICE directory detected at:\n  %s\n', juiceKernelDir);
    fprintf('If you download ESA JUICE mission kernels there, this script can\n');
    fprintf('also compare reconstructed states against the spacecraft SPK.\n');
end

fprintf('\nDone.\n');

function localLoadJuiceKernels(kernelDir)
    files = dir(fullfile(kernelDir, '**', '*'));
    for k = 1:numel(files)
        if files(k).isdir
            continue
        end
        [~,~,ext] = fileparts(files(k).name);
        if any(strcmpi(ext, {'.bsp','.bc','.tf','.tls','.tpc','.tsc','.tm','.bpc'}))
            try
                cspice_furnsh(fullfile(files(k).folder, files(k).name));
            catch
            end
        end
    end
end