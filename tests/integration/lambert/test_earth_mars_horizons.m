clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/astroToolbox')
startup

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
% Primary solver: Izzo
% -------------------------------------------------------------------------
solI = astro.lambert.solveIzzo(r1, r2, tof, sun.mu, false);

if ~solI.converged
    error('FAIL: Izzo solver did not converge.');
end

% -------------------------------------------------------------------------
% Optional comparison solver: Universal
% -------------------------------------------------------------------------
solU = astro.lambert.solveUniversal(r1, r2, tof, sun.mu, false);

if ~solU.converged
    warning('Universal solver did not converge. Continuing with Izzo-only validation.');
    dv1diff = NaN;
    dv2diff = NaN;
else
    dv1diff = norm(solU.v1 - solI.v1);
    dv2diff = norm(solU.v2 - solI.v2);
end

% -------------------------------------------------------------------------
% Propagation-based validation of Izzo solution
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

fprintf('Earth-Mars Horizons Lambert validation\n');
fprintf('--------------------------------------\n');
fprintf('Departure UTC              : %s\n', depUTC);
fprintf('Arrival UTC                : %s\n', arrUTC);
fprintf('Time of flight             : %.6f days\n', tof/86400);
fprintf('Izzo iterations            : %d\n', solI.iterations);
fprintf('Final position error       : %.6e km\n', posError);
fprintf('Final velocity error       : %.6e km/s\n', velError);

if ~isnan(dv1diff)
    fprintf('||v1_U - v1_I||            : %.6e km/s\n', dv1diff);
    fprintf('||v2_U - v2_I||            : %.6e km/s\n', dv2diff);
end

% -------------------------------------------------------------------------
% Tolerances
% -------------------------------------------------------------------------
posTol = 1e-3;   % km   = 1 meter
velTol = 1e-8;   % km/s

if posError < posTol && velError < velTol
    disp('PASS: Earth-Mars Izzo Lambert solution validated.')
else
    error('FAIL: Earth-Mars Izzo Lambert validation did not satisfy tolerances.')
end