clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/AstroToolbox')
startup

earth = astro.bodies.getBody('earth');

% Orbit
altitude = 700;
a = earth.radius + altitude;
e = 0;
inc = deg2rad(51.6);
RAAN = deg2rad(30);
omega = 0;

Torb = astro.maneuvers.orbitalPeriod(a, earth.mu);
t = linspace(0, Torb, 1200).';

% Station
lat = deg2rad(25.0);
lon0 = deg2rad(55.0);
hStation = 0.0;
minElevationDeg = 10.0;

omegaEarth = 7.2921159e-5;
theta0 = 0.0;

azDeg = zeros(numel(t),1);
elDeg = zeros(numel(t),1);
rngKm = zeros(numel(t),1);
visible = false(numel(t),1);

for k = 1:numel(t)
    M = 2*pi*t(k)/Torb;
    theta = M;

    [rSc, ~] = astro.coords.coe2rv(a, e, inc, RAAN, omega, theta, earth.mu);

    rStation = astro.visibility.stationECI( ...
        lat, lon0, hStation, earth.radius, t(k), theta0, omegaEarth);

    lonNow = lon0 + theta0 + omegaEarth*t(k);

    aer = astro.visibility.azimuthElevationRange(rSc, rStation, lat, lonNow);
    azDeg(k) = aer.azimuthDeg;
    elDeg(k) = aer.elevationDeg;
    rngKm(k) = aer.range;

    visible(k) = elDeg(k) >= minElevationDeg;
end

figure('Color','w');
plot(t/60, azDeg, 'LineWidth', 1.2);
grid on
xlabel('Time [min]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Azimuth [deg]', 'FontSize', 13, 'FontWeight', 'bold');
title('Topocentric Azimuth vs Time', 'FontSize', 15, 'FontWeight', 'bold');

figure('Color','w');
hold on
plot(t/60, elDeg, 'LineWidth', 1.2);
yline(minElevationDeg, '--', 'LineWidth', 1.2);
plot(t(visible)/60, elDeg(visible), '.', 'MarkerSize', 8);
grid on
xlabel('Time [min]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Elevation [deg]', 'FontSize', 13, 'FontWeight', 'bold');
title('Topocentric Elevation vs Time', 'FontSize', 15, 'FontWeight', 'bold');
legend('Elevation', 'Mask', 'Visible samples', 'Location', 'best');

figure('Color','w');
plot(t/60, rngKm, 'LineWidth', 1.2);
grid on
xlabel('Time [min]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Range [km]', 'FontSize', 13, 'FontWeight', 'bold');
title('Topocentric Range vs Time', 'FontSize', 15, 'FontWeight', 'bold');