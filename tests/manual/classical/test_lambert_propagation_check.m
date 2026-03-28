clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/astroToolbox')
startup

earth = astro.bodies.earth();

% ==========================================================
% Example: Lambert solution validated by two-body propagation
% ==========================================================

r1  = [5000; 10000; 2100];    % km
r2  = [-14600; 2500; 7000];   % km
tof = 3600;                   % s

% Solve Lambert
sol = astro.lambert.solveUniversal(r1, r2, tof, earth.mu, false);

fprintf('Lambert converged     : %d\n', sol.converged);
fprintf('Lambert iterations    : %d\n', sol.iterations);
fprintf('Lambert TOF error     : %.6e s\n', sol.tofError);

% Propagate from r1 with Lambert departure velocity
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

% Plot propagated trajectory and endpoint geometry
figure
plot3(out.x(:,1), out.x(:,2), out.x(:,3), 'LineWidth', 1.5)
hold on
plot3(r1(1), r1(2), r1(3), 'o', 'MarkerSize', 8, 'LineWidth', 1.5)
plot3(r2(1), r2(2), r2(3), 's', 'MarkerSize', 8, 'LineWidth', 1.5)
plot3(rf(1), rf(2), rf(3), 'x', 'MarkerSize', 10, 'LineWidth', 1.5)

grid on
axis equal
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
title('Lambert propagation validation')
legend('Propagated transfer arc', 'Initial position r_1', ...
       'Target position r_2', 'Propagated final position', ...
       'Location', 'best')