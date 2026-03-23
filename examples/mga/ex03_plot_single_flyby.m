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

rpVenus = venus.radius + 300;  % km

astro.plot.plotFlyby2D( ...
    traj.legs(1).vInfArr, ...
    traj.legs(2).vInfDep, ...
    venus.mu, rpVenus, +1, venus.radius, 'Venus');