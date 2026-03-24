clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/astroToolbox')
startup

% ==========================================================
% Earth -> Mars patched-conics mission design example
% Real heliocentric states from JPL Horizons
% Primary Lambert solver: Izzo
% ==========================================================

sun   = astro.bodies.getBody('sun');
earth = astro.bodies.getBody('earth');
mars  = astro.bodies.getBody('mars');

% ----------------------------------------------------------
% User inputs
% ----------------------------------------------------------
depUTC = '2026-11-12 00:00:00';
arrUTC = '2027-08-22 00:00:00';

hLEO = 300;   % km parking orbit altitude at Earth
hLMO = 300;   % km capture orbit altitude at Mars

rLEO = earth.radius + hLEO;
rLMO = mars.radius  + hLMO;

% Mean heliocentric orbital radii for SOI estimate
aEarth = 149597870.7;   % km
aMars  = 227939200.0;   % km

% ----------------------------------------------------------
% Fetch real heliocentric states
% ----------------------------------------------------------
earthState = astro.ephem.getHorizonsState('399', depUTC, '500@10');
marsState  = astro.ephem.getHorizonsState('499', arrUTC, '500@10');

r1 = earthState.r;
r2 = marsState.r;

t1 = datetime(depUTC, 'InputFormat', 'yyyy-MM-dd HH:mm:ss', 'TimeZone', 'UTC');
t2 = datetime(arrUTC, 'InputFormat', 'yyyy-MM-dd HH:mm:ss', 'TimeZone', 'UTC');
tof = seconds(t2 - t1);

% ----------------------------------------------------------
% Lambert transfer
% ----------------------------------------------------------
sol = astro.lambert.solveIzzo(r1, r2, tof, sun.mu, false);

if ~sol.converged
    error('Izzo solver did not converge for Earth-Mars patched-conics example.');
end

% ----------------------------------------------------------
% Hyperbolic excess velocities
% ----------------------------------------------------------
vInfDepVec = astro.maneuvers.hyperbolicExcess(sol.v1, earthState.v);
vInfArrVec = astro.maneuvers.hyperbolicExcess(sol.v2, marsState.v);

vInfDep = norm(vInfDepVec);
vInfArr = norm(vInfArrVec);

% ----------------------------------------------------------
% C3 and patched-conic departure/capture
% ----------------------------------------------------------
c3Dep = astro.maneuvers.C3(vInfDep);

dep = astro.maneuvers.departureDeltaV(earth.mu, rLEO, vInfDep);
arr = astro.maneuvers.arrivalDeltaV(mars.mu,  rLMO, vInfArr);

% ----------------------------------------------------------
% Sphere of influence estimates
% ----------------------------------------------------------
rSOI_Earth = astro.maneuvers.sphereOfInfluence(aEarth, earth.mu, sun.mu);
rSOI_Mars  = astro.maneuvers.sphereOfInfluence(aMars,  mars.mu,  sun.mu);

% ----------------------------------------------------------
% Propagation check
% ----------------------------------------------------------
x0 = [r1; sol.v1];

opts.RelTol = 1e-12;
opts.AbsTol = 1e-12;
opts.Solver = 'ode113';

out = astro.propagators.propagate( ...
    @(t,x) astro.propagators.eomTwoBody(t, x, sun.mu), ...
    [0 tof], x0, opts);

rf = out.x(end,1:3).';
vf = out.x(end,4:6).';

posError = norm(rf - r2);
velError = norm(vf - sol.v2);

% ----------------------------------------------------------
% Report
% ----------------------------------------------------------
fprintf('Earth-Mars patched-conics mission example\n');
fprintf('-----------------------------------------\n');
fprintf('Departure UTC                  : %s\n', depUTC);
fprintf('Arrival UTC                    : %s\n', arrUTC);
fprintf('Time of flight                 : %.6f days\n', tof/86400);
fprintf('\n');

fprintf('Lambert solver                 : Izzo\n');
fprintf('Converged                      : %d\n', sol.converged);
fprintf('Iterations                     : %d\n', sol.iterations);
fprintf('TOF residual                   : %.6e\n', sol.tofError);
fprintf('Propagation position error     : %.6e km\n', posError);
fprintf('Propagation velocity error     : %.6e km/s\n', velError);
fprintf('\n');

fprintf('Departure v_inf                : %.6f km/s\n', vInfDep);
fprintf('Arrival   v_inf                : %.6f km/s\n', vInfArr);
fprintf('Departure C3                   : %.6f km^2/s^2\n', c3Dep);
fprintf('\n');

fprintf('Earth parking orbit radius     : %.6f km\n', rLEO);
fprintf('Mars  capture orbit radius     : %.6f km\n', rLMO);
fprintf('Earth departure delta-v        : %.6f km/s\n', dep.deltaV);
fprintf('Mars  arrival delta-v          : %.6f km/s\n', arr.deltaV);
fprintf('\n');

fprintf('Earth SOI radius               : %.6f km\n', rSOI_Earth);
fprintf('Mars  SOI radius               : %.6f km\n', rSOI_Mars);

% ----------------------------------------------------------
% Plot heliocentric transfer arc
% ----------------------------------------------------------
figure
hold on
astro.plot.plotOrbit2D(out.x(:,1:3), 'LineWidth', 1.5);
plot(r1(1), r1(2), 'o', 'MarkerSize', 8, 'LineWidth', 1.5);
plot(r2(1), r2(2), 's', 'MarkerSize', 8, 'LineWidth', 1.5);
plot(rf(1), rf(2), 'x', 'MarkerSize', 10, 'LineWidth', 1.5);
plot(0, 0, 'o', 'MarkerSize', 10, 'LineWidth', 1.5);

xlabel('x [km]')
ylabel('y [km]')
title('Earth-to-Mars patched-conics transfer (heliocentric)')
legend('Lambert arc', 'Earth departure state', 'Mars arrival state', ...
       'Propagated final point', 'Sun', 'Location', 'best')
axis equal
grid on