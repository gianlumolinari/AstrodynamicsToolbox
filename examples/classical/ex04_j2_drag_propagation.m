clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/astroToolbox')
startup

earth = astro.bodies.getBody('earth');

% Initial orbit
a = earth.radius + 500;    % km
e = 0.001;
inc = deg2rad(51.6);
RAAN = deg2rad(40);
omega = deg2rad(20);
theta = deg2rad(0);

[r0, v0] = astro.coords.coe2rv(a, e, inc, RAAN, omega, theta, earth.mu);
x0 = [r0; v0];

% Perturbation settings
pert.useJ2 = true;
pert.useDrag = true;
pert.Cd = 2.2;
pert.AoverM = 0.01;              % m^2/kg
pert.rho0 = 3.614e-13;           % kg/m^3 at ~700 km-ish reference style placeholder
pert.H = 88.667;                 % km
pert.omegaBody = 7.2921159e-5;   % rad/s

Torb = astro.maneuvers.orbitalPeriod(a, earth.mu);
tspan = [0, 10*Torb];

opts.RelTol = 1e-11;
opts.AbsTol = 1e-11;
opts.Solver = 'ode113';

out = astro.propagators.propagate( ...
    @(t,x) astro.propagators.eomCowell(t, x, earth, pert), ...
    tspan, x0, opts);

figure('Color','w');
hold on
axis equal
grid on

ang = linspace(0, 2*pi, 300);
plot(earth.radius*cos(ang), earth.radius*sin(ang), 'k', 'LineWidth', 1.2);
plot(out.x(:,1), out.x(:,2), 'LineWidth', 1.2);

xlabel('x [km]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('y [km]', 'FontSize', 13, 'FontWeight', 'bold');
title('J2 + Drag Propagation Example', 'FontSize', 15, 'FontWeight', 'bold');
legend('Earth', 'Trajectory', 'Location', 'best');

fprintf('Propagation complete over %.2f orbits.\n', tspan(2)/Torb);