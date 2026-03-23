% Clear environment
clear; clc; close all;
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultAxesFontSize', 20);
set(0, 'DefaultLineLineWidth', 2);

% Load SPICE kernels
cspice_furnsh('/Users/gianlucamolinari/Documents/MATLAB/mice/Kernels/naif0012.tls');
cspice_furnsh('/Users/gianlucamolinari/Documents/MATLAB/mice/Kernels/de440s.bsp');
cspice_furnsh('/Users/gianlucamolinari/Documents/MATLAB/mice/Kernels/pck00010.tpc');

% Time span (1 year, 1-hour step)
startUTC = '2025-06-01T00:00:00';
endUTC   = '2026-06-01T00:00:00';

et_start = cspice_str2et(startUTC);
et_end   = cspice_str2et(endUTC);
dt = 3600;  % 1 hour in seconds

et_vector = et_start:dt:et_end;
N = length(et_vector);

% Preallocate position arrays
moon_eci     = zeros(3, N);  % Moon in ECI frame (w.r.t. Earth)
moon_synodic = zeros(3, N);  % Moon in rotating synodic frame (w.r.t. Earth)

for i = 1:N
    et = et_vector(i);

    % Get Moon and Earth positions in J2000 inertial frame (w.r.t. Sun)
    [moon_pos_sun,  ~] = cspice_spkpos('MOON',  et, 'J2000', 'NONE', 'SUN');
    [earth_pos_sun, ~] = cspice_spkpos('EARTH', et, 'J2000', 'NONE', 'SUN');

    % Relative Moon position in ECI (Moon - Earth)
    rel_moon = moon_pos_sun - earth_pos_sun;
    moon_eci(:, i) = rel_moon;

    % Build DCM for rotating Sun–Earth synodic frame
    x_hat = earth_pos_sun / norm(earth_pos_sun);  % from Sun to Earth
    x_hat = x_hat(:);

    earth_state = cspice_spkezr('EARTH', et, 'J2000', 'NONE', 'SUN');
    v_earth = earth_state(4:6);

    z_hat = cross(earth_pos_sun, v_earth);
    z_hat = z_hat / norm(z_hat);
    z_hat = z_hat(:);

    y_hat = cross(z_hat, x_hat);
    y_hat = y_hat / norm(y_hat);
    y_hat = y_hat(:);

    R = [x_hat, y_hat, z_hat];  % Rotation matrix: inertial → rotating

    % Transform Moon position into rotating frame
    moon_synodic(:, i) = R' * rel_moon;
end

%% Plot 1: Moon orbit in ECI frame

figure('Name', 'Moon Orbit in ECI Frame', 'NumberTitle', 'off');
plot3(moon_eci(1,:), moon_eci(2,:), moon_eci(3,:), 'r', 'DisplayName', 'Moon Orbit (ECI)');
hold on;
plot3(0, 0, 0, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 10, 'DisplayName', 'Earth');
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
title('Moon Orbit in J2000 Inertial Frame (ECI)');
grid on;
axis equal;
view(3);
legend('Location', 'best');

%% Plot 2: Moon orbit in rotating Sun–Earth synodic frame

figure('Name', 'Moon Orbit in Synodic Frame', 'NumberTitle', 'off');
plot3(moon_synodic(1,:), moon_synodic(2,:), moon_synodic(3,:), 'b', 'DisplayName', 'Moon Orbit (Synodic)');
hold on;
plot3(0, 0, 0, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 10, 'DisplayName', 'Earth (Origin)');
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
title('Moon Orbit in Sun Earth Synodic Frame');
grid on;
axis equal;
view(3);
legend('Location', 'best');

% Unload SPICE kernels
cspice_kclear;