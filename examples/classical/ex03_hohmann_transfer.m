clc;
clear;
close all;

startup

% ==========================================================
% Example: Hohmann transfer between two circular Earth orbits
% ==========================================================

earth = astro.bodies.earth();

r1 = 7000;    % km
r2 = 12000;   % km

man = astro.maneuvers.hohmann(earth.mu, r1, r2);

fprintf('Hohmann transfer from %.1f km to %.1f km\n', r1, r2);
fprintf('-----------------------------------------\n');
fprintf('Transfer semi-major axis  = %.6f km\n', man.aTrans);
fprintf('Initial circular speed    = %.6f km/s\n', man.vCirc1);
fprintf('Final circular speed      = %.6f km/s\n', man.vCirc2);
fprintf('Transfer periapsis speed  = %.6f km/s\n', man.vPeriTrans);
fprintf('Transfer apoapsis speed   = %.6f km/s\n', man.vApoTrans);
fprintf('First impulse dv1         = %.6f km/s\n', man.dv1);
fprintf('Second impulse dv2        = %.6f km/s\n', man.dv2);
fprintf('Total delta-v             = %.6f km/s\n', man.dvTot);
fprintf('Transfer time of flight   = %.6f s\n', man.tof);
fprintf('Transfer time of flight   = %.6f hr\n', man.tof/3600);

% Simple 2D geometry plot
theta = linspace(0, 2*pi, 400);

x1 = r1 * cos(theta);
y1 = r1 * sin(theta);

x2 = r2 * cos(theta);
y2 = r2 * sin(theta);

aT = man.aTrans;
eT = abs(r2 - r1) / (r1 + r2);
pT = aT * (1 - eT^2);

thetaT = linspace(0, pi, 300);
rT = pT ./ (1 + eT*cos(thetaT));

xT = rT .* cos(thetaT);
yT = rT .* sin(thetaT);

figure
plot(0, 0, 'o', 'MarkerSize', 8, 'LineWidth', 1.5)
hold on
plot(x1, y1, 'LineWidth', 1.5)
plot(x2, y2, 'LineWidth', 1.5)
plot(xT, yT, 'LineWidth', 1.8)
grid on
axis equal
xlabel('x [km]')
ylabel('y [km]')
title('Hohmann transfer between circular orbits')
legend('Earth center', 'Initial orbit', 'Final orbit', 'Transfer arc', ...
    'Location', 'best')