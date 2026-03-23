clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/AstroToolbox')
startup
astro.ephem.loadSpiceKernels(fullfile(pwd, 'data', 'spice'));

sun = astro.bodies.getBody('sun');

sequence = {'earth','venus','earth','jupiter'};

dates = {
    '2026-01-01 00:00:00'
    '2026-06-01 00:00:00'
    '2027-01-01 00:00:00'
    '2029-01-01 00:00:00'
};

traj = astro.mga.propagateMGA(sequence, dates, sun.mu, 'spice');

bodyMap.venus = astro.bodies.getBody('venus');
bodyMap.earth = astro.bodies.getBody('earth');

minAltMap.venus = 300;   % km
minAltMap.earth = 300;   % km

report = astro.mga.evaluateMGA(traj, bodyMap, minAltMap, 1e-3);

disp('MGA leg summary:')
for k = 1:numel(traj.legs)
    fprintf('\nLeg %d: %s -> %s\n', k, upper(traj.legs(k).bodyDepart), upper(traj.legs(k).bodyArrive));
    fprintf('  TOF                 : %.2f days\n', traj.legs(k).tofDays);
    fprintf('  |v_inf,dep|         : %.6f km/s\n', traj.legs(k).vInfDepMag);
    fprintf('  |v_inf,arr|         : %.6f km/s\n', traj.legs(k).vInfArrMag);
end

fprintf('\nFlyby feasibility summary\n');
fprintf('-------------------------\n');

for k = 1:report.nFlybys
    fb = report.flybys(k);
    chk = fb.check;

    fprintf('\nFlyby at %s\n', upper(fb.body));
    fprintf('  Incoming |v_inf^-|       : %.6f km/s\n', chk.magIn);
    fprintf('  Required |v_inf^+|       : %.6f km/s\n', chk.magOutReq);
    fprintf('  Magnitude mismatch       : %.6f km/s\n', chk.magMismatch);
    fprintf('  Required turn angle      : %.6f deg\n', chk.deltaReqDeg);
    fprintf('  Max turn angle @ rpMin   : %.6f deg\n', chk.deltaMaxDeg);
    fprintf('  Min altitude             : %.2f km\n', fb.minAltitude);
    fprintf('  Feasible                 : %d\n', chk.feasible);
end