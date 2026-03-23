clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/AstroToolbox')
startup
astro.ephem.loadSpiceKernels(fullfile(pwd, 'data', 'spice'));

sun = astro.bodies.getBody('sun');
venus = astro.bodies.getBody('venus');
earth = astro.bodies.getBody('earth');

sequence = {'earth','venus','earth','jupiter'};

dates = {
    '2026-01-01 00:00:00'
    '2026-06-01 00:00:00'
    '2027-01-01 00:00:00'
    '2029-01-01 00:00:00'
};

traj = astro.mga.propagateMGA(sequence, dates, sun.mu, 'spice');

disp('MGA summary:')
for k = 1:numel(traj.legs)
    fprintf('\nLeg %d: %s -> %s\n', k, upper(traj.legs(k).bodyDepart), upper(traj.legs(k).bodyArrive));
    fprintf('  TOF           : %.2f days\n', traj.legs(k).tofDays);
    fprintf('  |v_inf,dep|   : %.6f km/s\n', traj.legs(k).vInfDepMag);
    fprintf('  |v_inf,arr|   : %.6f km/s\n', traj.legs(k).vInfArrMag);
end

% Example planar flyby turn at Venus using the arrival v_infinity from leg 1
rpVenus = venus.radius + 300;  % km
flybyV = astro.mga.flybyTurn(traj.legs(1).vInfArr, venus.mu, rpVenus, +1);

fprintf('\nExample Venus flyby turn\n');
fprintf('  rp            : %.2f km\n', rpVenus);
fprintf('  |v_inf|       : %.6f km/s\n', flybyV.vInfMag);
fprintf('  turn angle    : %.6f deg\n', flybyV.deltaDeg);

% Compare outgoing flyby v_inf to required departure v_inf of next leg
vInfReqNext = traj.legs(2).vInfDep;
fprintf('  Required next-leg |v_inf,dep| : %.6f km/s\n', norm(vInfReqNext));
fprintf('  Example flyby out |v_inf,out| : %.6f km/s\n', norm(flybyV.vInfOut));
fprintf('  Vector mismatch magnitude      : %.6f km/s\n', norm(flybyV.vInfOut - vInfReqNext));