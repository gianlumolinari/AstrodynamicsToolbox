clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/astroToolbox')
startup

earth = astro.bodies.earth();

r1 = [5000; 10000; 2100];   % km
r2 = [-14600; 2500; 7000];  % km
tof = 3600;                 % s

sol = astro.lambert.solveUniversal(r1, r2, tof, earth.mu, false);

fprintf('Converged  : %d\n', sol.converged);
fprintf('Iterations : %d\n', sol.iterations);
fprintf('TOF error  : %.6e s\n', sol.tofError);

if sol.converged && abs(sol.tofError) < 1e-8
    disp('PASS: Lambert solver converged.')
else
    error('FAIL: Lambert solver did not converge satisfactorily.')
end