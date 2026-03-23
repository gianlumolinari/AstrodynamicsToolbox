clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/AstroToolbox')
startup
astro.ephem.loadSpiceKernels(fullfile(pwd, 'data', 'spice'));

sun = astro.bodies.getBody('sun');

% ----------------------------------------------------------
% MGA sequence
% ----------------------------------------------------------
sequence = {'earth','venus','earth','jupiter'};
t0UTC = '2026-01-01 00:00:00';

% Initial guess [days]
x0 = [260; 110; 1200];

% Bounds [days]
lb = [80; 80; 300];
ub = [400; 400; 1600];

% Flyby body data
bodyMap.venus = astro.bodies.getBody('venus');
bodyMap.earth = astro.bodies.getBody('earth');

minAltMap.venus = 300;   % km
minAltMap.earth = 300;   % km

obj = @(x) astro.mga.mgaObjective(x, sequence, t0UTC, sun.mu, bodyMap, minAltMap, 'spice');
nonlcon = @(x) astro.mga.mgaConstraints(x, sequence, t0UTC, sun.mu, bodyMap, minAltMap, 'spice');

fprintf('Running constrained MGA optimization...\n');

opts = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'Algorithm', 'sqp', ...
    'MaxFunctionEvaluations', 4000, ...
    'StepTolerance', 1e-8, ...
    'ConstraintTolerance', 1e-4);

[xOpt, fOpt, exitflag, output] = fmincon(obj, x0, [], [], [], [], lb, ub, nonlcon, opts);

fprintf('\nOptimization complete.\n');
fprintf('Exit flag: %d\n', exitflag);
fprintf('Objective : %.6f\n', fOpt);

datesOpt = astro.mga.unpackDatesFromTOF(t0UTC, xOpt);
trajOpt = astro.mga.propagateMGA(sequence, datesOpt, sun.mu, 'spice');
reportOpt = astro.mga.evaluateMGA(trajOpt, bodyMap, minAltMap, 1e-3);

fprintf('\nOptimized encounter dates:\n');
for k = 1:numel(datesOpt)
    fprintf('  %s : %s\n', upper(sequence{k}), datesOpt{k});
end

fprintf('\nOptimized flyby feasibility summary\n');
fprintf('-----------------------------------\n');
for k = 1:reportOpt.nFlybys
    fb = reportOpt.flybys(k);
    chk = fb.check;

    fprintf('\nFlyby at %s\n', upper(fb.body));
    fprintf('  Incoming |v_inf^-|       : %.6f km/s\n', chk.magIn);
    fprintf('  Required |v_inf^+|       : %.6f km/s\n', chk.magOutReq);
    fprintf('  Magnitude mismatch       : %.6f km/s\n', chk.magMismatch);
    fprintf('  Required turn angle      : %.6f deg\n', chk.deltaReqDeg);
    fprintf('  Max turn angle @ rpMin   : %.6f deg\n', chk.deltaMaxDeg);
    fprintf('  Feasible                 : %d\n', chk.feasible);
end

% ----------------------------------------------------------
% Flyby plot at Venus
% ----------------------------------------------------------
venus = astro.bodies.getBody('venus');
rpVenus = venus.radius + minAltMap.venus;

astro.plot.plotFlyby2D( ...
    trajOpt.legs(1).vInfArr, ...
    trajOpt.legs(2).vInfDep, ...
    venus.mu, rpVenus, +1, venus.radius, 'Venus');

% ----------------------------------------------------------
% Heliocentric overview
% ----------------------------------------------------------
figure('Color','w');
hold on
axis equal
grid on

plot(0, 0, 'o', 'MarkerSize', 10, 'LineWidth', 1.5);

for k = 1:numel(trajOpt.legs)
    x0Leg = [trajOpt.legs(k).r1; trajOpt.legs(k).v1];
    tspan = [0, trajOpt.legs(k).tofSec];

    optsProp.RelTol = 1e-11;
    optsProp.AbsTol = 1e-11;
    optsProp.Solver = 'ode113';

    out = astro.propagators.propagate( ...
        @(t,x) astro.propagators.eomTwoBody(t, x, sun.mu), ...
        tspan, x0Leg, optsProp);

    astro.plot.plotOrbit2D(out.x(:,1:3), 'LineWidth', 1.3);
    plot(trajOpt.legs(k).r1(1), trajOpt.legs(k).r1(2), 'o', 'MarkerSize', 7, 'LineWidth', 1.2);
    plot(trajOpt.legs(k).r2(1), trajOpt.legs(k).r2(2), 's', 'MarkerSize', 7, 'LineWidth', 1.2);
end

xlabel('x [km]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('y [km]', 'FontSize', 13, 'FontWeight', 'bold');
title('Constrained MGA Heliocentric Transfer Overview', ...
    'FontSize', 15, 'FontWeight', 'bold');
legend('Sun', 'Leg arcs / encounters', 'Location', 'best');