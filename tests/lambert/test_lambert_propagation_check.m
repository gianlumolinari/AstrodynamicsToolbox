clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/astroToolbox')
startup

earth = astro.bodies.earth();

r1  = [5000; 10000; 2100];    % km
r2  = [-14600; 2500; 7000];   % km
tof = 3600;                   % s

sol = astro.lambert.solveUniversal(r1, r2, tof, earth.mu, false);

if ~sol.converged
    error('FAIL: Lambert solver did not converge.')
end

x0 = [r1; sol.v1];

opts.RelTol = 1e-12;
opts.AbsTol = 1e-12;
opts.Solver = 'ode113';

out = astro.propagators.propagate( ...
    @(t,x) astro.propagators.eomTwoBody(t, x, earth.mu), ...
    [0 tof], x0, opts);

rf = out.x(end,1:3).';
vf = out.x(end,4:6).';

posError = norm(rf - r2);
velError = norm(vf - sol.v2);

fprintf('Final position error  = %.6e km\n', posError);
fprintf('Final velocity error  = %.6e km/s\n', velError);

posTol = 1e-6;    % km
velTol = 1e-9;    % km/s

if posError < posTol && velError < velTol
    disp('PASS: Lambert solution validated by propagation.')
else
    error('FAIL: Lambert propagation validation did not satisfy tolerances.')
end