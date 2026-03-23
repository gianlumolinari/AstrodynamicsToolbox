clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/astroToolbox')
startup
astro.ephem.loadSpiceKernels(fullfile(pwd, 'data', 'spice'));

% ==========================================================
% Example: Earth eclipse check using a simple cylindrical shadow model
% Plot in the Earth-centered frame for clarity
% ==========================================================

sun   = astro.bodies.getBody('sun');
earth = astro.bodies.getBody('earth');

epochUTC = '2026-11-12 00:00:00';

% Earth heliocentric state from SPICE
earthState = astro.ephem.getState('spice', 'earth', epochUTC);

rEarth = earthState.r;
rSun = [0; 0; 0];   % Sun-centered frame

% Spacecraft circular orbit around Earth
altitude = 700;                          % km
rOrbit = earth.radius + altitude;        % km

theta = linspace(0, 2*pi, 400);
inShadow = false(size(theta));

% Store orbit in Earth-centered frame
xRel = zeros(size(theta));
yRel = zeros(size(theta));

for k = 1:numel(theta)
    % Earth-centered spacecraft position
    rRel = rOrbit * [cos(theta(k)); sin(theta(k)); 0];
    rSc = rEarth + rRel;   % heliocentric position for eclipse check
    
    out = astro.geometry.isInShadow(rSc, rSun, rEarth, earth.radius);
    inShadow(k) = out.inShadow;
    
    xRel(k) = rRel(1);
    yRel(k) = rRel(2);
end

% ----------------------------------------------------------
% Compute Sun direction in Earth-centered frame
% ----------------------------------------------------------
sunDirFromEarth = rSun - rEarth;
sunHat = sunDirFromEarth / norm(sunDirFromEarth);
shadowHat = -sunHat;

shadowLength = 3 * rOrbit;
shadowLine = [zeros(3,1), shadowLength * shadowHat];

% ----------------------------------------------------------
% Plot in Earth-centered frame
% ----------------------------------------------------------
figure('Color','w');
hold on
axis equal
grid on

% Draw Earth as a filled circle
ang = linspace(0, 2*pi, 300);
xEarth = earth.radius * cos(ang);
yEarth = earth.radius * sin(ang);
fill(xEarth, yEarth, [0.85 0.92 1.0], 'EdgeColor', 'k', 'LineWidth', 1.2);

% Plot spacecraft orbit arcs
plot(xRel(~inShadow), yRel(~inShadow), 'LineWidth', 1.5);
plot(xRel(inShadow),  yRel(inShadow),  '.', 'MarkerSize', 10);

% Plot Sun direction arrow
quiver(0, 0, sunHat(1)*1.6*rOrbit, sunHat(2)*1.6*rOrbit, 0, ...
    'LineWidth', 1.5, 'MaxHeadSize', 0.5);

% Plot shadow axis
plot(shadowLine(1,:), shadowLine(2,:), '--', 'LineWidth', 1.5);

xlabel('x_{ECI-like} [km]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('y_{ECI-like} [km]', 'FontSize', 13, 'FontWeight', 'bold');
title('Earth Orbit Eclipse Check (Cylindrical Shadow Model)', ...
    'FontSize', 15, 'FontWeight', 'bold');

legend('Earth', 'Sunlit arc', 'Eclipsed arc', 'Sun direction', 'Shadow axis', ...
    'Location', 'best');

% Axis limits
lim = 1.8 * rOrbit;
xlim([-lim, lim]);
ylim([-lim, lim]);

fprintf('Epoch: %s\n', epochUTC);
fprintf('Orbit altitude: %.2f km\n', altitude);
fprintf('Fraction of sampled orbit in eclipse: %.4f\n', mean(inShadow));