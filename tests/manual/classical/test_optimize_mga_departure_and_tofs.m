clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/AstroToolbox')
startup
astro.ephem.loadSpiceKernels(fullfile(pwd, 'data', 'spice'));

sun = astro.bodies.getBody('sun');

sequence = {'earth','venus','earth','jupiter'};
t0RefUTC = '2026-01-01 00:00:00';

% x = [departureOffsetDays; tof1; tof2; tof3]
x0 = [0; 260; 110; 1200];

lb = [-400; 80; 80; 300];
ub = [ 400; 500; 500; 1800];

bodyMap.venus = astro.bodies.getBody('venus');
bodyMap.earth = astro.bodies.getBody('earth');

minAltMap.venus = 300;
minAltMap.earth = 300;

obj = @(x) localObjective(x, sequence, t0RefUTC, sun.mu, bodyMap, minAltMap);
nonlcon = @(x) localConstraints(x, sequence, t0RefUTC, sun.mu, bodyMap, minAltMap);

opts = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'Algorithm', 'sqp', ...
    'MaxFunctionEvaluations', 6000, ...
    'StepTolerance', 1e-8, ...
    'ConstraintTolerance', 1e-4);

fprintf('Running constrained MGA optimization with variable departure date...\n');
[xOpt, fOpt, exitflag] = fmincon(obj, x0, [], [], [], [], lb, ub, nonlcon, opts);

fprintf('\nOptimization complete.\n');
fprintf('Exit flag: %d\n', exitflag);
fprintf('Objective : %.6f\n', fOpt);

datesOpt = astro.mga.unpackDatesFromDecisionVector(t0RefUTC, xOpt);
trajOpt = astro.mga.propagateMGA(sequence, datesOpt, sun.mu, 'spice');
reportOpt = astro.mga.evaluateMGA(trajOpt, bodyMap, minAltMap, 1e-3);

fprintf('\nOptimized encounter dates:\n');
for k = 1:numel(datesOpt)
    fprintf('  %s : %s\n', upper(sequence{k}), datesOpt{k});
end

fprintf('\nFlyby feasibility summary\n');
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

function J = localObjective(x, sequence, t0RefUTC, muSun, bodyMap, minAltMap)
    dates = astro.mga.unpackDatesFromDecisionVector(t0RefUTC, x);
    tofDays = x(2:end);
    J = astro.mga.mgaObjective(tofDays, sequence, dates{1}, muSun, bodyMap, minAltMap, 'spice');
end

function [c, ceq] = localConstraints(x, sequence, t0RefUTC, muSun, bodyMap, minAltMap)
    dates = astro.mga.unpackDatesFromDecisionVector(t0RefUTC, x);
    tofDays = x(2:end);
    [c, ceq] = astro.mga.mgaConstraints(tofDays, sequence, dates{1}, muSun, bodyMap, minAltMap, 'spice');
end