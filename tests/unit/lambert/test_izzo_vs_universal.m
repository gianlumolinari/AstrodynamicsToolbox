clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/astroToolbox')
startup

earth = astro.bodies.earth();

r1  = [5000; 10000; 2100];    % km
r2  = [-14600; 2500; 7000];   % km
tof = 3600;                   % s

solU = astro.lambert.solveUniversal(r1, r2, tof, earth.mu, false);
solI = astro.lambert.solveIzzo(r1, r2, tof, earth.mu, false);

if ~solU.converged
    error('FAIL: Universal solver did not converge.')
end

if ~solI.converged
    error('FAIL: Izzo solver did not converge.')
end

dv1diff = norm(solU.v1 - solI.v1);
dv2diff = norm(solU.v2 - solI.v2);

fprintf('Difference in departure velocity = %.6e km/s\n', dv1diff);
fprintf('Difference in arrival velocity   = %.6e km/s\n', dv2diff);
fprintf('Izzo iterations                  = %d\n', solI.iterations);

tol = 1e-10;

if dv1diff < tol && dv2diff < tol
    disp('PASS: Izzo solver matches universal solver within tolerance.')
else
    error('FAIL: Izzo solver does not match universal solver within tolerance.')
end