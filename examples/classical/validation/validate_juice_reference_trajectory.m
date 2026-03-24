clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/AstroToolbox')
startup

% Load generic planetary kernels
astro.ephem.loadSpiceKernels(fullfile(pwd, 'data', 'spice'));

% ==========================================================
% JUICE Reference Trajectory Reconstruction & SPICE Truth Comparison
% High-fidelity propagated reconstruction
% ==========================================================
sun     = astro.bodies.getBody('sun');
earth   = astro.bodies.getBody('earth');
moon    = astro.bodies.getBody('moon');
venus   = astro.bodies.getBody('venus');
mars    = astro.bodies.getBody('mars');
jupiter = astro.bodies.getBody('jupiter');
saturn  = astro.bodies.getBody('saturn');

% ----------------------------------------------------------
% Heliocentric sequence & Dates
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
planetNames = {'Earth','Earth (LEGA)','Venus','Earth','Earth','Jupiter'};

fprintf('\n============================================================\n');
fprintf('JUICE Reference Trajectory Reconstruction & Analysis\n');
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
% Flyby feasibility setup & Leg Summary
% ----------------------------------------------------------
bodyMap.earth = earth;
bodyMap.venus = venus;
minAltMap.earth = 300;   % km
minAltMap.venus = 300;   % km
report = astro.mga.evaluateMGA(traj, bodyMap, minAltMap, 1e-3);

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

launchC3 = traj.legs(1).vInfDepMag^2;
arrivalVInfJupiter = traj.legs(end).vInfArrMag;
totalTOFdays = sum([traj.legs.tofDays]);

fprintf('\nMission-level summary:\n');
fprintf('  Launch C3                 : %.6f km^2/s^2\n', launchC3);
fprintf('  Jupiter arrival v_inf     : %.6f km/s\n', arrivalVInfJupiter);
fprintf('  Total time of flight      : %.2f days\n', totalTOFdays);
fprintf('  Total time of flight      : %.2f years\n', totalTOFdays / 365.25);

% ----------------------------------------------------------
% Flyby 2D Plots
% ----------------------------------------------------------
rpEarth = earth.radius + minAltMap.earth;
rpVenus = venus.radius + minAltMap.venus;

astro.plot.plotFlyby2D(traj.legs(1).vInfArr, traj.legs(2).vInfDep, earth.mu, rpEarth, +1, earth.radius, 'Earth (LEGA)');
astro.plot.plotFlyby2D(traj.legs(2).vInfArr, traj.legs(3).vInfDep, venus.mu, rpVenus, +1, venus.radius, 'Venus');
astro.plot.plotFlyby2D(traj.legs(3).vInfArr, traj.legs(4).vInfDep, earth.mu, rpEarth, +1, earth.radius, 'Earth');
astro.plot.plotFlyby2D(traj.legs(4).vInfArr, traj.legs(5).vInfDep, earth.mu, rpEarth, +1, earth.radius, 'Earth');

% ----------------------------------------------------------
% Build high-fidelity propagated leg samples for plotting
% ----------------------------------------------------------
configHF = struct();
configHF.muSun = sun.mu;
configHF.bodyNames = { ...
    'VENUS BARYCENTER', ...
    'EARTH', ...
    'MOON', ...
    'MARS BARYCENTER', ...
    'JUPITER BARYCENTER', ...
    'SATURN BARYCENTER'};

configHF.muBodies = [ ...
    venus.mu, ...
    earth.mu, ...
    moon.mu, ...
    mars.mu, ...
    jupiter.mu, ...
    saturn.mu];

% JUICE SRP parameters are uncertain here; keep easy to toggle
configHF.useSRP   = false;
configHF.Cr       = 1.3;
configHF.A_over_m = 0.01;

configHF.RelTol = 1e-11;
configHF.AbsTol = 1e-11;

rReconHF = [];
rReconLegStart = zeros(numel(traj.legs),3);
rReconLegEnd   = zeros(numel(traj.legs),3);

fprintf('\n============================================================\n');
fprintf('High-Fidelity Leg Propagation\n');
fprintf('============================================================\n');

for k = 1:numel(traj.legs)
    x0Leg = [traj.legs(k).r1; traj.legs(k).v1];
    tfLeg = traj.legs(k).tofSec;
    et0Leg = cspice_str2et(dates{k});

    configLeg = configHF;
    configLeg.et0 = et0Leg;

    outLeg = astro.propagators.propagateHighFidelity(x0Leg, et0Leg, tfLeg, configLeg);

    rReconHF = [rReconHF; outLeg.x(:,1:3)]; %#ok<AGROW>
    rReconLegStart(k,:) = outLeg.x(1,1:3);
    rReconLegEnd(k,:)   = outLeg.x(end,1:3);

    fprintf('Leg %d propagated: %s -> %s, TOF = %.2f days\n', ...
        k, upper(traj.legs(k).bodyDepart), upper(traj.legs(k).bodyArrive), traj.legs(k).tofDays);
end

% ----------------------------------------------------------
% SPICE Truth Comparison & Heliocentric Overview Plot
% ----------------------------------------------------------
juiceKernelDir = fullfile(pwd, 'data', 'spice', 'juice');
hasJuiceDir = exist(juiceKernelDir, 'dir') == 7;

figure('Color','w');
hold on
axis equal
grid on

hSun   = plot(0, 0, 'ko', 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', 'Sun');
hRecon = plot(rReconHF(:,1), rReconHF(:,2), 'b-', 'LineWidth', 1.3, ...
    'DisplayName', 'High-fidelity propagated reconstruction');

xlabel('x [km]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('y [km]', 'FontSize', 13, 'FontWeight', 'bold');

if hasJuiceDir
    fprintf('\n============================================================\n');
    fprintf('JUICE SPICE Truth Validation\n');
    fprintf('============================================================\n');
    fprintf('Loading JUICE mission kernels from:\n  %s\n', juiceKernelDir);

    try
        localLoadJuiceKernels(juiceKernelDir);

        candidateTargets = {'JUICE', 'JUICE SPACECRAFT', '-28'};
        juiceTarget = '';
        testEpoch = dates{1};

        for k = 1:numel(candidateTargets)
            try
                tmp = astro.ephem.getSpiceState(candidateTargets{k}, testEpoch, 'SUN', 'J2000', 'NONE');
                juiceTarget = tmp.target;
                break
            catch
            end
        end

        if isempty(juiceTarget)
            error('Could not resolve JUICE target name automatically.');
        end

        fprintf('Resolved JUICE spacecraft target as: %s\n', juiceTarget);

        t0 = datetime(dates{1}, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
        tf = datetime(dates{end}, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
        nSamples = 1200;

        tHist = linspace(posixtime(t0), posixtime(tf), nSamples);
        tHist = datetime(tHist, 'ConvertFrom', 'posixtime');

        rJuice = zeros(nSamples, 3);

        fprintf('Sampling real JUICE trajectory from SPICE...\n');
        for k = 1:nSamples
            epochUTC = datestr(tHist(k), 'yyyy-mm-dd HH:MM:SS');
            sc = astro.ephem.getSpiceState(juiceTarget, epochUTC, 'SUN', 'J2000', 'NONE');
            rJuice(k,:) = sc.r(:).';
        end

        hTruth    = plot(rJuice(:,1), rJuice(:,2), 'g--', 'LineWidth', 1.4, 'DisplayName', 'JUICE SPICE truth');
        hTruthPts = plot(nan, nan, 'go', 'MarkerSize', 7, 'LineWidth', 1.2, 'DisplayName', 'JUICE encounter states');
        hReconPts = plot(nan, nan, 'bs', 'MarkerSize', 7, 'LineWidth', 1.2, 'DisplayName', 'Reconstructed encounter states');

        fprintf('\nEncounter state comparison:\n');
        fprintf('------------------------------------------------------------\n');
        for k = 1:numel(dates)
            epochUTC = dates{k};
            sc = astro.ephem.getSpiceState(juiceTarget, epochUTC, 'SUN', 'J2000', 'NONE');
            plot(sc.r(1), sc.r(2), 'go', 'MarkerSize', 7, 'LineWidth', 1.2);

            if k == 1
                rRef = traj.legs(1).r1;
                vRef = traj.legs(1).v1;
            elseif k == numel(dates)
                rRef = traj.legs(end).r2;
                vRef = traj.legs(end).v2;
            else
                rRef = traj.legs(k-1).r2;
                vRef = traj.legs(k-1).v2;
            end

            plot(rRef(1), rRef(2), 'bs', 'MarkerSize', 7, 'LineWidth', 1.2);

            text(sc.r(1), sc.r(2), ['  ' planetNames{k}], ...
                'FontSize', 10, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');

            posErr = norm(sc.r - rRef);
            velErr = norm(sc.v - vRef);

            fprintf('%-19s  pos err = %12.3f km   vel err = %10.6f km/s\n', ...
                epochUTC, posErr, velErr);
        end

        title('JUICE: SPICE Truth vs High-Fidelity Reconstruction', ...
            'FontSize', 15, 'FontWeight', 'bold');
        legend([hSun, hRecon, hTruth, hTruthPts, hReconPts], 'Location', 'best');

    catch ME
        fprintf('\nJUICE validation block failed.\nReason: %s\n', ME.message);
        title('JUICE Reference Trajectory: Heliocentric Overview', ...
            'FontSize', 15, 'FontWeight', 'bold');
        legend([hSun, hRecon], 'Location', 'best');
    end
else
    fprintf('\nNo JUICE mission SPICE directory detected at:\n  %s\n', juiceKernelDir);
    fprintf('If you download ESA JUICE mission kernels there, this script will\n');
    fprintf('automatically compare reconstructed states against the spacecraft SPK.\n');

    hReconPts = plot(nan, nan, 'bs', 'MarkerSize', 7, 'LineWidth', 1.2, ...
        'DisplayName', 'Encounter states');

    for k = 1:numel(traj.legs)
        plot(traj.legs(k).r1(1), traj.legs(k).r1(2), 'bs', 'MarkerSize', 7, 'LineWidth', 1.2);
        plot(traj.legs(k).r2(1), traj.legs(k).r2(2), 'bs', 'MarkerSize', 7, 'LineWidth', 1.2);
    end

    title('JUICE Reference Trajectory: High-Fidelity Heliocentric Overview', ...
        'FontSize', 15, 'FontWeight', 'bold');
    legend([hSun, hRecon, hReconPts], 'Location', 'best');
end

fprintf('\nDone.\n');

% ==========================================================
% Helper Functions
% ==========================================================
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