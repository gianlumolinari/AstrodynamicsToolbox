clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/astroToolbox')
startup

% ==========================================================
% Example: simple ground-station access for a circular orbit
% Earth-centered inertial-style demonstration
% ==========================================================

earth = astro.bodies.getBody('earth');

% Orbit
altitude = 700;                       % km
a = earth.radius + altitude;          % km
e = 0;
inc = deg2rad(51.6);
RAAN = deg2rad(30);
omega = 0;

Torb = astro.maneuvers.orbitalPeriod(a, earth.mu);
t = linspace(0, Torb, 1000).';

% Ground station (spherical Earth approximation, fixed in inertial frame
% for this simple first example)
lat = deg2rad(25.0);    % deg
lon = deg2rad(55.0);    % deg
rStation = earth.radius * [cos(lat)*cos(lon); cos(lat)*sin(lon); sin(lat)];

% Spacecraft states
rHist = zeros(numel(t), 3);
elevDeg = zeros(numel(t), 1);
visible = false(numel(t), 1);

minElevationDeg = 10;

for k = 1:numel(t)
    M = 2*pi*t(k)/Torb;
    theta = M;  % circular orbit
    
    [rSc, ~] = astro.coords.coe2rv(a, e, inc, RAAN, omega, theta, earth.mu);
    rHist(k,:) = rSc(:).';
    
    vis = astro.visibility.isVisibleFromStation(rSc, rStation, minElevationDeg);
    elevDeg(k) = vis.elevation.elevationDeg;
    visible(k) = vis.isVisible;
end

intervals = astro.visibility.accessIntervals(t, visible);

% ----------------------------------------------------------
% Plot 1: Orbit and station direction
% ----------------------------------------------------------
figure('Color','w');
hold on
grid on
axis equal

ang = linspace(0, 2*pi, 300);
plot(earth.radius*cos(ang), earth.radius*sin(ang), 'k', 'LineWidth', 1.2);
plot(rHist(:,1), rHist(:,2), 'LineWidth', 1.2);
plot(rHist(visible,1), rHist(visible,2), '.', 'MarkerSize', 8);
plot(rStation(1), rStation(2), 'o', 'MarkerSize', 8, 'LineWidth', 1.5);

xlabel('x [km]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('y [km]', 'FontSize', 13, 'FontWeight', 'bold');
title('Ground Station Access Example', 'FontSize', 15, 'FontWeight', 'bold');
legend('Earth', 'Orbit', 'Visible arc', 'Ground station', 'Location', 'best');

% ----------------------------------------------------------
% Plot 2: Elevation vs time
% ----------------------------------------------------------
figure('Color','w');
hold on
grid on

plot(t/60, elevDeg, 'LineWidth', 1.3);
yline(minElevationDeg, '--', 'LineWidth', 1.2);

xlabel('Time [min]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Elevation [deg]', 'FontSize', 13, 'FontWeight', 'bold');
title('Elevation Angle vs Time', 'FontSize', 15, 'FontWeight', 'bold');
legend('Elevation', 'Mask angle', 'Location', 'best');

fprintf('Orbit period: %.2f min\n', Torb/60);
fprintf('Elevation mask: %.2f deg\n', minElevationDeg);
fprintf('Number of access intervals: %d\n', numel(intervals));

for k = 1:numel(intervals)
    durSec = intervals(k).duration;
    fprintf('  Interval %d: start = %.2f min, stop = %.2f min, duration = %.2f min\n', ...
        k, intervals(k).start/60, intervals(k).stop/60, durSec/60);
end