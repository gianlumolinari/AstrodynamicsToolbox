clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/astroToolbox')
startup

% ==========================================================
% Earth -> Mars Lambert transfer using JPL Horizons states
% Primary solver: Izzo
% Universal solver kept only as a consistency check
% ==========================================================

sun = astro.bodies.getBody('sun');

depUTC = '2026-11-12 00:00:00';
arrUTC = '2027-08-22 00:00:00';

earthState = astro.ephem.getHorizonsState('399', depUTC, '500@10');
marsState  = astro.ephem.getHorizonsState('499', arrUTC, '500@10');

r1 = earthState.r;
r2 = marsState.r;

t1 = datetime(depUTC, 'InputFormat', 'yyyy-MM-dd HH:mm:ss', 'TimeZone', 'UTC');
t2 = datetime(arrUTC, 'InputFormat', 'yyyy-MM-dd HH:mm:ss', 'TimeZone', 'UTC');
tof = seconds(t2 - t1);

% -------------------------------------------------------------------------
% Primary Lambert solution: Izzo
% -------------------------------------------------------------------------
solI = astro.lambert.solveIzzo(r1, r2, tof, sun.mu, false);

if ~solI.converged
    error('Izzo solver did not converge for the Earth-Mars case.');
end

% -------------------------------------------------------------------------
% Optional cross-check: Universal solver
% -------------------------------------------------------------------------
solU = astro.lambert.solveUniversal(r1, r2, tof, sun.mu, false);

if solU.converged
    dv1diff = norm(solU.v1 - solI.v1);
    dv2diff = norm(solU.v2 - solI.v2);
else
    dv1diff = NaN;
    dv2diff = NaN;
end

% -------------------------------------------------------------------------
% Propagation-based validation using Izzo result
% -------------------------------------------------------------------------
x0 = [r1; solI.v1];

opts.RelTol = 1e-12;
opts.AbsTol = 1e-12;
opts.Solver = 'ode113';

out = astro.propagators.propagate( ...
    @(t,x) astro.propagators.eomTwoBody(t, x, sun.mu), ...
    [0 tof], x0, opts);

rf = out.x(end,1:3).';
vf = out.x(end,4:6).';

posError = norm(rf - r2);
velError = norm(vf - solI.v2);

% -------------------------------------------------------------------------
% Mission-style quantities
% -------------------------------------------------------------------------
vInfDep = norm(solI.v1 - earthState.v);
vInfArr = norm(solI.v2 - marsState.v);

fprintf('Earth-Mars Lambert transfer using JPL Horizons states\n');
fprintf('----------------------------------------------------\n');
fprintf('Departure UTC                 : %s\n', depUTC);
fprintf('Arrival UTC                   : %s\n', arrUTC);
fprintf('Time of flight                : %.6f days\n', tof/86400);
fprintf('\n');

fprintf('Izzo converged                : %d\n', solI.converged);
fprintf('Izzo iterations               : %d\n', solI.iterations);
fprintf('Izzo TOF error                : %.6e\n', solI.tofError);

if ~isnan(dv1diff)
    fprintf('||v1_U - v1_I||               : %.6e km/s\n', dv1diff);
    fprintf('||v2_U - v2_I||               : %.6e km/s\n', dv2diff);
else
    fprintf('Universal solver              : did not converge (ignored)\n');
end

fprintf('\n');
fprintf('Propagated final position err : %.6e km\n', posError);
fprintf('Propagated final velocity err : %.6e km/s\n', velError);
fprintf('Departure v_inf wrt Earth     : %.6f km/s\n', vInfDep);
fprintf('Arrival   v_inf wrt Mars      : %.6f km/s\n', vInfArr);

% -------------------------------------------------------------------------
% Plot transfer arc in heliocentric ecliptic plane
% -------------------------------------------------------------------------
figure;
hold on;

astro.plot.plotOrbit2D(out.x(:,1:3), 'LineWidth', 1.5);
plot(r1(1), r1(2), 'o', 'MarkerSize', 8, 'LineWidth', 1.5);
plot(r2(1), r2(2), 's', 'MarkerSize', 8, 'LineWidth', 1.5);
plot(rf(1), rf(2), 'x', 'MarkerSize', 10, 'LineWidth', 1.5);
plot(0, 0, 'o', 'MarkerSize', 10, 'LineWidth', 1.5);

xlabel('x [km]');
ylabel('y [km]');
title('Earth-to-Mars Lambert transfer (JPL Horizons states)');
legend('Lambert arc', 'Earth at departure', 'Mars at arrival', ...
       'Propagated final point', 'Sun', 'Location', 'best');
axis equal;
grid on;