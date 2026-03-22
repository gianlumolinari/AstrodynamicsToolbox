clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/astroToolbox')
startup

earth = astro.bodies.earth();

% Test geometry
r1  = [5000; 10000; 2100];    % km
r2  = [-14600; 2500; 7000];   % km
tof = 3600;                   % s

% Solve with both interfaces
solU = astro.lambert.solveUniversal(r1, r2, tof, earth.mu, false);
solI = astro.lambert.solveIzzo(r1, r2, tof, earth.mu, false);

fprintf('=== Universal Solver ===\n');
fprintf('Converged   : %d\n', solU.converged);
fprintf('Iterations  : %d\n', solU.iterations);
fprintf('TOF error   : %.6e s\n\n', solU.tofError);

fprintf('=== Izzo Interface ===\n');
fprintf('Converged   : %d\n', solI.converged);
fprintf('Iterations  : %d\n', solI.iterations);
fprintf('TOF error   : %.6e s\n', solI.tofError);
fprintf('Backend     : %s\n\n', solI.backend);

dv1diff = norm(solU.v1 - solI.v1);
dv2diff = norm(solU.v2 - solI.v2);

fprintf('||v1_U - v1_I|| = %.6e km/s\n', dv1diff);
fprintf('||v2_U - v2_I|| = %.6e km/s\n', dv2diff);

% Propagation check using Izzo interface output
x0 = [r1; solI.v1];

opts.RelTol = 1e-12;
opts.AbsTol = 1e-12;
opts.Solver = 'ode113';

out = astro.propagators.propagate( ...
    @(t,x) astro.propagators.eomTwoBody(t, x, earth.mu), ...
    [0 tof], x0, opts);

rf = out.x(end,1:3).';
vf = out.x(end,4:6).';

posError = norm(rf - r2);
velError = norm(vf - solI.v2);

fprintf('Izzo-interface propagated position error = %.6e km\n', posError);
fprintf('Izzo-interface propagated velocity error = %.6e km/s\n', velError);

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
title('Universal vs Izzo-interface Lambert comparison')
legend('Propagated arc', 'r_1', 'r_2', 'Propagated final point', ...
    'Location', 'best')