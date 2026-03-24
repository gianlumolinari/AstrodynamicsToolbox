clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/astroToolbox')
startup
astro.ephem.loadSpiceKernels(fullfile(pwd, 'data', 'spice'));

sun   = astro.bodies.getBody('sun');
earth = astro.bodies.getBody('earth');

epochUTC = '2026-11-12 00:00:00';
earthState = astro.ephem.getState('spice', 'earth', epochUTC);

rEarth = earthState.r;
rSun = [0;0;0];

altitude = 700;
rOrbit = earth.radius + altitude;

theta = linspace(0, 2*pi, 500);
stateCode = strings(size(theta));
xRel = zeros(size(theta));
yRel = zeros(size(theta));

for k = 1:numel(theta)
    rRel = rOrbit * [cos(theta(k)); sin(theta(k)); 0];
    rSc = rEarth + rRel;

    out = astro.geometry.isInShadowConical(rSc, rSun, rEarth, earth.radius, sun.radius);
    stateCode(k) = string(out.state);

    xRel(k) = rRel(1);
    yRel(k) = rRel(2);
end

isUmbra = stateCode == "umbra";
isPen = stateCode == "penumbra";
isSun = stateCode == "sunlight";

sunDirFromEarth = rSun - rEarth;
sunHat = sunDirFromEarth / norm(sunDirFromEarth);
shadowHat = -sunHat;

figure('Color','w');
hold on
axis equal
grid on

ang = linspace(0, 2*pi, 300);
fill(earth.radius*cos(ang), earth.radius*sin(ang), [0.85 0.92 1.0], ...
    'EdgeColor', 'k', 'LineWidth', 1.2);

plot(xRel(isSun), yRel(isSun), '.', 'MarkerSize', 8);
plot(xRel(isPen), yRel(isPen), '.', 'MarkerSize', 8);
plot(xRel(isUmbra), yRel(isUmbra), '.', 'MarkerSize', 8);

quiver(0,0,sunHat(1)*1.6*rOrbit,sunHat(2)*1.6*rOrbit,0,'LineWidth',1.5);
quiver(0,0,shadowHat(1)*1.6*rOrbit,shadowHat(2)*1.6*rOrbit,0,'LineWidth',1.5,'LineStyle','--');

xlabel('x [km]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('y [km]', 'FontSize', 13, 'FontWeight', 'bold');
title('Conical Eclipse Model: Sunlight / Penumbra / Umbra', ...
    'FontSize', 15, 'FontWeight', 'bold');
legend('Earth', 'Sunlit arc', 'Penumbra arc', 'Umbra arc', ...
    'Sun direction', 'Shadow axis', 'Location', 'best');

fprintf('Epoch: %s\n', epochUTC);
fprintf('Sunlit fraction   : %.4f\n', mean(isSun));
fprintf('Penumbra fraction : %.4f\n', mean(isPen));
fprintf('Umbra fraction    : %.4f\n', mean(isUmbra));