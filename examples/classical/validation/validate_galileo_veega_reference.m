clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/AstroToolbox')
startup
astro.ephem.loadSpiceKernels(fullfile(pwd, 'data', 'spice'));

% ==========================================================
% Galileo VEEGA reference reconstruction
%
% Verified historical sequence:
%   Earth launch   : 1989-10-18
%   Venus flyby    : 1990-02-10
%   Earth flyby 1  : 1990-12-08
%   Earth flyby 2  : 1992-12-08
%   Jupiter arrival: 1995-12-07
%
% This is a patched-conics / Lambert reconstruction intended
% as a polished toolbox capability demonstration, not an
% exact recreation of the flown navigation solution.
% ==========================================================

sun     = astro.bodies.getBody('sun');
venus   = astro.bodies.getBody('venus');
earth   = astro.bodies.getBody('earth');
jupiter = astro.bodies.getBody('jupiter');

% ----------------------------------------------------------
% Historical Galileo VEEGA sequence
% ----------------------------------------------------------
sequence = {'earth','venus','earth','earth','jupiter'};

dates = {
    '1989-10-18 00:00:00'
    '1990-02-10 00:00:00'
    '1990-12-08 00:00:00'
    '1992-12-08 00:00:00'
    '1995-12-07 00:00:00'
};

% ----------------------------------------------------------
% Propagate MGA legs
% ----------------------------------------------------------
traj = astro.mga.propagateMGA(sequence, dates, sun.mu, 'spice');

bodyMap.venus = venus;
bodyMap.earth = earth;

minAltMap.venus = 300;   % km
minAltMap.earth = 300;   % km

report = astro.mga.evaluateMGA(traj, bodyMap, minAltMap, 1e-3);

% ----------------------------------------------------------
% Summary
% ----------------------------------------------------------
fprintf('\n');
fprintf('============================================================\n');
fprintf('Galileo VEEGA Reference Reconstruction\n');
fprintf('============================================================\n');

fprintf('\nEncounter dates used:\n');
for k = 1:numel(sequence)
    fprintf('  %-8s : %s\n', upper(sequence{k}), dates{k});
end

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

legHandles = gobjects(numel(traj.legs),1);
depHandles = gobjects(numel(traj.legs),1);
arrHandles = gobjects(numel(traj.legs),1);

for k = 1:numel(traj.legs)
    x0Leg = [traj.legs(k).r1; traj.legs(k).v1];
    tspan = [0, traj.legs(k).tofSec];

    optsProp.RelTol = 1e-11;
    optsProp.AbsTol = 1e-11;
    optsProp.Solver = 'ode113';

    out = astro.propagators.propagate( ...
        @(t,x) astro.propagators.eomTwoBody(t, x, sun.mu), ...
        tspan, x0Leg, optsProp);

    legHandles(k) = plot(out.x(:,1), out.x(:,2), 'LineWidth', 1.4);

    depHandles(k) = plot(traj.legs(k).r1(1), traj.legs(k).r1(2), ...
        'o', 'MarkerSize', 7, 'LineWidth', 1.2);

    arrHandles(k) = plot(traj.legs(k).r2(1), traj.legs(k).r2(2), ...
        's', 'MarkerSize', 7, 'LineWidth', 1.2);
end

xlabel('x [km]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('y [km]', 'FontSize', 13, 'FontWeight', 'bold');
title('Galileo VEEGA Reference Reconstruction: Heliocentric Overview', ...
    'FontSize', 15, 'FontWeight', 'bold');
legend([plot(nan,nan,'o','MarkerSize',10,'LineWidth',1.5), ...
        legHandles(1), depHandles(1), arrHandles(1)], ...
       {'Sun','Lambert leg','Departure encounter','Arrival encounter'}, ...
       'Location','best');

% ----------------------------------------------------------
% Flyby geometry plots
% ----------------------------------------------------------
rpVenus = venus.radius + minAltMap.venus;
rpEarth = earth.radius + minAltMap.earth;

% Venus flyby: leg 1 arrival vs leg 2 departure
astro.plot.plotFlyby2D( ...
    traj.legs(1).vInfArr, ...
    traj.legs(2).vInfDep, ...
    venus.mu, rpVenus, +1, venus.radius, 'Venus');

% Earth flyby 1: leg 2 arrival vs leg 3 departure
astro.plot.plotFlyby2D( ...
    traj.legs(2).vInfArr, ...
    traj.legs(3).vInfDep, ...
    earth.mu, rpEarth, +1, earth.radius, 'Earth');

% Earth flyby 2: leg 3 arrival vs leg 4 departure
astro.plot.plotFlyby2D( ...
    traj.legs(3).vInfArr, ...
    traj.legs(4).vInfDep, ...
    earth.mu, rpEarth, +1, earth.radius, 'Earth');

% ----------------------------------------------------------
% Mission-level scalar summary
% ----------------------------------------------------------
launchC3 = traj.legs(1).vInfDepMag^2;
arrivalVInfJupiter = traj.legs(end).vInfArrMag;
totalTOFdays = sum([traj.legs.tofDays]);

fprintf('\nMission-level summary:\n');
fprintf('  Launch C3                 : %.6f km^2/s^2\n', launchC3);
fprintf('  Jupiter arrival v_inf     : %.6f km/s\n', arrivalVInfJupiter);
fprintf('  Total time of flight      : %.2f days\n', totalTOFdays);
fprintf('  Total time of flight      : %.2f years\n', totalTOFdays / 365.25);

fprintf('\nDone.\n');