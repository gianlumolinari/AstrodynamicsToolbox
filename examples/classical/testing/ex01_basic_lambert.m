clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/astroToolbox')
startup

earth = astro.bodies.earth();

% Example geometry
r1 = [5000; 10000; 2100];   % km
r2 = [-14600; 2500; 7000];  % km
tof = 3600;                 % s

sol = astro.lambert.solveUniversal(r1, r2, tof, earth.mu, false);

disp('=== Lambert Solution ===')
fprintf('Converged     : %d\n', sol.converged);
fprintf('Iterations    : %d\n', sol.iterations);
fprintf('z             : %.12f\n', sol.z);
fprintf('TOF error     : %.6e s\n', sol.tofError);

disp('Departure velocity v1 [km/s]:')
disp(sol.v1)

disp('Arrival velocity v2 [km/s]:')
disp(sol.v2)

% Simple geometry plot
figure
plot3(r1(1), r1(2), r1(3), 'o', 'MarkerSize', 8, 'LineWidth', 1.5)
hold on
plot3(r2(1), r2(2), r2(3), 's', 'MarkerSize', 8, 'LineWidth', 1.5)
plot3([0 r1(1)], [0 r1(2)], [0 r1(3)], '--', 'LineWidth', 1.0)
plot3([0 r2(1)], [0 r2(2)], [0 r2(3)], '--', 'LineWidth', 1.0)
grid on
axis equal
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
title('Lambert transfer geometry')
legend('r_1', 'r_2', 'Initial radius', 'Final radius', 'Location', 'best')