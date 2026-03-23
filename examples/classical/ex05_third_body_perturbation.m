clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/AstroToolbox')
startup
astro.ephem.loadSpiceKernels(fullfile(pwd, 'data', 'spice'));

earth = astro.bodies.getBody('earth');
moon  = astro.bodies.getBody('moon');

% ----------------------------------------------------------
% Spacecraft initial orbit about Earth
% ----------------------------------------------------------
a = earth.radius + 20000;   % km
e = 0.1;
inc = deg2rad(30);
RAAN = deg2rad(20);
omega = deg2rad(40);
theta = deg2rad(0);

[r0, v0] = astro.coords.coe2rv(a, e, inc, RAAN, omega, theta, earth.mu);
x0 = [r0; v0];

% ----------------------------------------------------------
% Moon state relative to Earth at reference epoch
% ----------------------------------------------------------
epochUTC = '2026-11-12 00:00:00';
moonState = astro.ephem.getSpiceState('MOON', epochUTC, 'EARTH', 'J2000', 'NONE');
rMoon = moonState.r;

% ----------------------------------------------------------
% Propagation with frozen third-body perturbation
% ----------------------------------------------------------
Torb = astro.maneuvers.orbitalPeriod(a, earth.mu);
tspan = [0, 5*Torb];

opts.RelTol = 1e-11;
opts.AbsTol = 1e-11;
opts.Solver = 'ode113';

rhs = @(t,x) localEOM(t, x, earth, moon, rMoon);
out = astro.propagators.propagate(rhs, tspan, x0, opts);

% ----------------------------------------------------------
% Plot 1: local Earth-centered zoom
% ----------------------------------------------------------
figure('Color','w');
hold on
axis equal
grid on

ang = linspace(0, 2*pi, 300);
plot(earth.radius*cos(ang), earth.radius*sin(ang), 'k', 'LineWidth', 1.2);
plot(out.x(:,1), out.x(:,2), 'LineWidth', 1.4);

xlabel('x [km]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('y [km]', 'FontSize', 13, 'FontWeight', 'bold');
title('Third-Body Perturbation Example: Earth-Centered Zoom', ...
    'FontSize', 15, 'FontWeight', 'bold');
legend('Earth', 'Trajectory', 'Location', 'best');

% ----------------------------------------------------------
% Plot 2: wide context including Moon
% ----------------------------------------------------------
figure('Color','w');
hold on
axis equal
grid on

plot(earth.radius*cos(ang), earth.radius*sin(ang), 'k', 'LineWidth', 1.2);
plot(out.x(:,1), out.x(:,2), 'LineWidth', 1.2);
plot(rMoon(1), rMoon(2), 'o', 'MarkerSize', 8, 'LineWidth', 1.5);

% Draw dashed line from Earth to Moon
plot([0 rMoon(1)], [0 rMoon(2)], '--', 'LineWidth', 1.0);

xlabel('x [km]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('y [km]', 'FontSize', 13, 'FontWeight', 'bold');
title('Third-Body Perturbation Example: Earth-Moon Context', ...
    'FontSize', 15, 'FontWeight', 'bold');
legend('Earth', 'Trajectory', 'Moon reference position', 'Earth-Moon line', ...
    'Location', 'best');

fprintf('Reference epoch: %s\n', epochUTC);
fprintf('Moon distance from Earth: %.2f km\n', norm(rMoon));
fprintf('Propagation complete over %.2f orbits.\n', tspan(2)/Torb);

function dx = localEOM(~, x, earth, moon, rMoon)
    r = x(1:3);
    v = x(4:6);

    rmag = norm(r);
    a2b = -earth.mu * r / rmag^3;
    a3b = astro.perturbations.accelThirdBody(r, rMoon, moon.mu);

    dx = [v; a2b + a3b];
end