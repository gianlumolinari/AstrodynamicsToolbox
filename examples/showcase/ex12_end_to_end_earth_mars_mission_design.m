clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/AstroToolbox')
startup
astro.ephem.loadSpiceKernels(fullfile(pwd, 'data', 'spice'));

% ==========================================================
% End-to-end Earth-Mars mission design showcase
%
% Demonstrates:
%   1) SPICE ephemerides
%   2) Lambert transfer design
%   3) Propagation-based Lambert validation
%   4) Patched-conics metrics (v_inf, C3, delta-v)
%   5) Heliocentric transfer plotting
%   6) Earth parking-orbit eclipse analysis
%   7) Ground-station visibility analysis
%   8) J2 + drag parking-orbit propagation
%   9) 3D central-body plotting utility
% ==========================================================

% ----------------------------------------------------------
% Bodies
% ----------------------------------------------------------
sun   = astro.bodies.getBody('sun');
earth = astro.bodies.getBody('earth');
mars  = astro.bodies.getBody('mars');

% ----------------------------------------------------------
% Mission dates
% ----------------------------------------------------------
depUTC = '2026-11-12 00:00:00';
arrUTC = '2027-08-22 00:00:00';

tDep = datetime(depUTC, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
tArr = datetime(arrUTC, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
tofSec = seconds(tArr - tDep);
tofDays = tofSec / 86400;

fprintf('\n');
fprintf('============================================================\n');
fprintf('End-to-End Earth-Mars Mission Design Showcase\n');
fprintf('============================================================\n');
fprintf('\nMission dates:\n');
fprintf('  Departure UTC : %s\n', depUTC);
fprintf('  Arrival UTC   : %s\n', arrUTC);
fprintf('  TOF           : %.2f days\n', tofDays);
fprintf('  TOF           : %.2f years\n', tofDays/365.25);

% ----------------------------------------------------------
% SPICE ephemerides
% ----------------------------------------------------------
earthState = astro.ephem.getState('spice', 'earth', depUTC);
marsState  = astro.ephem.getState('spice', 'mars',  arrUTC);

r1 = earthState.r;
vE = earthState.v;

r2 = marsState.r;
vM = marsState.v;

% ----------------------------------------------------------
% Lambert solve
% ----------------------------------------------------------
sol = astro.lambert.solveIzzo(r1, r2, tofSec, sun.mu, false);

if ~sol.converged
    error('Lambert solver did not converge.');
end

fprintf('\nLambert solution:\n');
fprintf('  Converged     : %d\n', sol.converged);
fprintf('  Iterations    : %d\n', sol.iterations);
fprintf('  TOF error     : %.6e s\n', sol.tofError);

v1 = sol.v1;
v2 = sol.v2;

% ----------------------------------------------------------
% Propagation-based validation
% ----------------------------------------------------------
x0 = [r1; v1];
optsProp.RelTol = 1e-12;
optsProp.AbsTol = 1e-12;
optsProp.Solver = 'ode113';

outLam = astro.propagators.propagate( ...
    @(t,x) astro.propagators.eomTwoBody(t, x, sun.mu), ...
    [0, tofSec], x0, optsProp);

rFinalProp = outLam.x(end,1:3).';
vFinalProp = outLam.x(end,4:6).';

posErr = norm(rFinalProp - r2);
velErr = norm(vFinalProp - v2);

fprintf('\nPropagation validation:\n');
fprintf('  Final position error : %.6e km\n', posErr);
fprintf('  Final velocity error : %.6e km/s\n', velErr);

% ----------------------------------------------------------
% Patched-conics mission metrics
% ----------------------------------------------------------
vInfDep = norm(v1 - vE);
vInfArr = norm(vM - v2);
C3 = vInfDep^2;

% Parking orbit assumptions
hLEO = 300;   % km
hLMO = 300;   % km
rLEO = earth.radius + hLEO;
rLMO = mars.radius + hLMO;

dep = astro.maneuvers.departureDeltaV(earth.mu, rLEO, vInfDep);
arr = astro.maneuvers.arrivalDeltaV(mars.mu, rLMO, vInfArr);

dvDep = dep.deltaV;
dvArr = arr.deltaV;
dvTot = dvDep + dvArr;

fprintf('\nPatched-conics summary:\n');
fprintf('  Departure v_inf       : %.6f km/s\n', vInfDep);
fprintf('  Arrival   v_inf       : %.6f km/s\n', vInfArr);
fprintf('  Launch C3             : %.6f km^2/s^2\n', C3);
fprintf('  Earth departure dv    : %.6f km/s\n', dvDep);
fprintf('  Mars  arrival dv      : %.6f km/s\n', dvArr);
fprintf('  Total mission dv      : %.6f km/s\n', dvTot);

% ----------------------------------------------------------
% Plot styling helpers
% ----------------------------------------------------------
ang = linspace(0, 2*pi, 400);
earthFaceColor = [0.78 0.86 0.96];
earthEdgeColor = [0.10 0.10 0.10];
orbitBlue = [0.00 0.20 0.95];
visibleGreen = [0.20 0.85 0.20];
penColor = [0.85 0.00 0.75];
umbraRed = [0.95 0.10 0.10];
stationMagenta = [0.85 0.00 0.85];

% ----------------------------------------------------------
% Heliocentric transfer plot
% ----------------------------------------------------------
figure('Color','w');
hold on
axis equal
grid on
box on

plot(0, 0, 'ko', 'MarkerSize', 10, 'LineWidth', 1.5);
plot(outLam.x(:,1), outLam.x(:,2), '-', 'Color', orbitBlue, 'LineWidth', 1.6);
plot(r1(1), r1(2), 'o', 'MarkerSize', 8, 'LineWidth', 1.4, ...
    'MarkerEdgeColor', visibleGreen, 'MarkerFaceColor', 'none');
plot(r2(1), r2(2), 's', 'MarkerSize', 8, 'LineWidth', 1.4, ...
    'MarkerEdgeColor', umbraRed, 'MarkerFaceColor', 'none');

text(r1(1), r1(2), '  Earth departure', 'FontSize', 10, 'FontWeight', 'bold');
text(r2(1), r2(2), '  Mars arrival', 'FontSize', 10, 'FontWeight', 'bold');

xlabel('x [km]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('y [km]', 'FontSize', 13, 'FontWeight', 'bold');
title('Earth-Mars Heliocentric Transfer', 'FontSize', 15, 'FontWeight', 'bold');
legend('Sun', 'Lambert transfer', 'Earth departure', 'Mars arrival', 'Location', 'best');

xAnn = min(outLam.x(:,1)) + 0.08*(max(outLam.x(:,1)) - min(outLam.x(:,1)));
yAnn = min(outLam.x(:,2)) + 0.15*(max(outLam.x(:,2)) - min(outLam.x(:,2)));
summaryText = sprintf(['TOF = %.1f days\n' ...
                       'C3 = %.3f km^2/s^2\n' ...
                       'v_{\\infty,dep} = %.3f km/s\n' ...
                       'v_{\\infty,arr} = %.3f km/s\n' ...
                       '\\Delta v_{tot} = %.3f km/s'], ...
                       tofDays, C3, vInfDep, vInfArr, dvTot);
text(xAnn, yAnn, summaryText, ...
    'FontSize', 10, ...
    'BackgroundColor', 'w', ...
    'EdgeColor', [0.3 0.3 0.3], ...
    'Margin', 8);

% ----------------------------------------------------------
% Earth parking-orbit eclipse analysis
% ----------------------------------------------------------
theta = linspace(0, 2*pi, 500);
stateCode = strings(size(theta));
xRel = zeros(size(theta));
yRel = zeros(size(theta));

rEarthHelio = earthState.r;
rSun = [0;0;0];

for k = 1:numel(theta)
    rRel = rLEO * [cos(theta(k)); sin(theta(k)); 0];
    rScHelio = rEarthHelio + rRel;

    sh = astro.geometry.isInShadowConical(rScHelio, rSun, rEarthHelio, earth.radius, sun.radius);
    stateCode(k) = string(sh.state);

    xRel(k) = rRel(1);
    yRel(k) = rRel(2);
end

isUmbra = stateCode == "umbra";
isPen   = stateCode == "penumbra";
isSun   = stateCode == "sunlight";

sunDirFromEarth = rSun - rEarthHelio;
sunHat = sunDirFromEarth / norm(sunDirFromEarth);
shadowHat = -sunHat;

figure('Color','w');
hold on
axis equal
grid on
box on

fill(earth.radius*cos(ang), earth.radius*sin(ang), earthFaceColor, ...
    'EdgeColor', earthEdgeColor, 'LineWidth', 1.3);

plot(xRel(isSun), yRel(isSun), '.', 'Color', orbitBlue, 'MarkerSize', 8);
plot(xRel(isPen), yRel(isPen), '.', 'Color', penColor, 'MarkerSize', 8);
plot(xRel(isUmbra), yRel(isUmbra), '.', 'Color', umbraRed, 'MarkerSize', 8);

quiver(0,0,sunHat(1)*1.6*rLEO,sunHat(2)*1.6*rLEO,0, ...
    'Color', 'k', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
quiver(0,0,shadowHat(1)*1.6*rLEO,shadowHat(2)*1.6*rLEO,0, ...
    'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);

xlabel('x [km]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('y [km]', 'FontSize', 13, 'FontWeight', 'bold');
title('Earth Parking Orbit Eclipse Classification', 'FontSize', 15, 'FontWeight', 'bold');
legend('Earth', 'Sunlit arc', 'Penumbra arc', 'Umbra arc', ...
    'Sun direction', 'Shadow axis', 'Location', 'best');

fprintf('\nEarth parking-orbit eclipse summary:\n');
fprintf('  Sunlit fraction   : %.4f\n', mean(isSun));
fprintf('  Penumbra fraction : %.4f\n', mean(isPen));
fprintf('  Umbra fraction    : %.4f\n', mean(isUmbra));

% ----------------------------------------------------------
% Ground-station visibility on parking orbit
% ----------------------------------------------------------
lat = deg2rad(25.0);      % Abu Dhabi-like
lon0 = deg2rad(55.0);
hStation = 0.0;
minElevationDeg = 10.0;

omegaEarth = 7.2921159e-5;   % rad/s
theta0 = 0.0;

Tpark = astro.maneuvers.orbitalPeriod(rLEO, earth.mu);
tVis = linspace(0, Tpark, 1200).';

elevDeg = zeros(numel(tVis),1);
visible = false(numel(tVis),1);
rScVisHist = zeros(numel(tVis),3);
rStationHist = zeros(numel(tVis),3);

for k = 1:numel(tVis)
    M = 2*pi*tVis(k)/Tpark;
    thetaNow = M;

    [rSc, ~] = astro.coords.coe2rv(rLEO, 0, deg2rad(28.5), 0, 0, thetaNow, earth.mu);
    rScVisHist(k,:) = rSc(:).';

    rStation = astro.visibility.stationECI( ...
        lat, lon0, hStation, earth.radius, tVis(k), theta0, omegaEarth);
    rStationHist(k,:) = rStation(:).';

    vis = astro.visibility.isVisibleFromStation(rSc, rStation, minElevationDeg);
    elevDeg(k) = vis.elevation.elevationDeg;
    visible(k) = vis.isVisible;
end

intervals = astro.visibility.accessIntervals(tVis, visible);

figure('Color','w');
hold on
grid on
axis equal
box on

plot(earth.radius*cos(ang), earth.radius*sin(ang), 'Color', earthEdgeColor, 'LineWidth', 1.2);
plot(rScVisHist(:,1), rScVisHist(:,2), '-', 'Color', orbitBlue, 'LineWidth', 1.2);
plot(rScVisHist(visible,1), rScVisHist(visible,2), '.', 'Color', visibleGreen, 'MarkerSize', 8);
plot(rStationHist(:,1), rStationHist(:,2), '--', 'Color', stationMagenta, 'LineWidth', 1.2);
plot(rStationHist(1,1), rStationHist(1,2), 'ko', 'MarkerSize', 8, 'LineWidth', 1.5);

xlabel('x [km]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('y [km]', 'FontSize', 13, 'FontWeight', 'bold');
title('Ground Station Access in Earth Parking Orbit', 'FontSize', 15, 'FontWeight', 'bold');
legend('Earth', 'Parking orbit', 'Visible arc', 'Station track in ECI', ...
    'Station at t_0', 'Location', 'best');

figure('Color','w');
hold on
grid on
box on

plot(tVis/60, elevDeg, '-', 'Color', orbitBlue, 'LineWidth', 1.3);
yline(minElevationDeg, '--', 'Color', umbraRed, 'LineWidth', 1.2);
plot(tVis(visible)/60, elevDeg(visible), '.', 'Color', visibleGreen, 'MarkerSize', 8);

xlabel('Time [min]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Elevation [deg]', 'FontSize', 13, 'FontWeight', 'bold');
title('Ground Station Elevation vs Time', 'FontSize', 15, 'FontWeight', 'bold');
legend('Elevation', 'Mask angle', 'Visible samples', 'Location', 'best');

% Optional 3D visibility view
figure('Color','w');
hold on
grid on
box on
axis equal

astro.plot.plotCentralBody(earth, earthFaceColor);
plot3(rScVisHist(:,1), rScVisHist(:,2), rScVisHist(:,3), ...
    '-', 'Color', orbitBlue, 'LineWidth', 1.2);
plot3(rScVisHist(visible,1), rScVisHist(visible,2), rScVisHist(visible,3), ...
    '.', 'Color', visibleGreen, 'MarkerSize', 8);
plot3(rStationHist(:,1), rStationHist(:,2), rStationHist(:,3), ...
    '--', 'Color', stationMagenta, 'LineWidth', 1.2);

xlabel('x [km]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('y [km]', 'FontSize', 13, 'FontWeight', 'bold');
zlabel('z [km]', 'FontSize', 13, 'FontWeight', 'bold');
title('Ground Station Access in Earth Parking Orbit (3D View)', ...
    'FontSize', 15, 'FontWeight', 'bold');
legend('Earth', 'Parking orbit', 'Visible arc', 'Station track in ECI', ...
    'Location', 'best');
view(35, 25)

fprintf('\nGround-station access summary:\n');
fprintf('  Elevation mask       : %.2f deg\n', minElevationDeg);
fprintf('  Number of intervals  : %d\n', numel(intervals));

for k = 1:numel(intervals)
    fprintf('  Interval %d: start = %.2f min, stop = %.2f min, duration = %.2f min\n', ...
        k, intervals(k).start/60, intervals(k).stop/60, intervals(k).duration/60);
end

% ----------------------------------------------------------
% J2 + drag propagation in parking orbit
% ----------------------------------------------------------
[r0Park, v0Park] = astro.coords.coe2rv(rLEO, 0.001, deg2rad(28.5), deg2rad(10), deg2rad(20), 0, earth.mu);
x0Park = [r0Park; v0Park];

pert.useJ2 = true;
pert.useDrag = true;
pert.Cd = 2.2;              % -
pert.AoverM = 0.01;         % m^2/kg
pert.rho0 = 3.614e-13;      % kg/m^3
pert.H = 88.667;            % km
pert.omegaBody = 7.2921159e-5;

tspanCowell = [0, 5*Tpark];
optsCowell.RelTol = 1e-11;
optsCowell.AbsTol = 1e-11;
optsCowell.Solver = 'ode113';

outCowell = astro.propagators.propagate( ...
    @(t,x) astro.propagators.eomCowell(t, x, earth, pert), ...
    tspanCowell, x0Park, optsCowell);

figure('Color','w');
hold on
grid on
box on
axis equal

astro.plot.plotCentralBody(earth, earthFaceColor);
plot3(outCowell.x(:,1), outCowell.x(:,2), outCowell.x(:,3), ...
    '-', 'Color', orbitBlue, 'LineWidth', 1.6);

xlabel('x [km]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('y [km]', 'FontSize', 13, 'FontWeight', 'bold');
zlabel('z [km]', 'FontSize', 13, 'FontWeight', 'bold');
title('Parking Orbit Propagation with J2 + Drag', 'FontSize', 15, 'FontWeight', 'bold');
legend('Earth', 'Trajectory', 'Location', 'best');
view(35, 25)

% ----------------------------------------------------------
% Final consolidated summary
% ----------------------------------------------------------
fprintf('\n');
fprintf('============================================================\n');
fprintf('Showcase summary\n');
fprintf('============================================================\n');
fprintf('Departure UTC             : %s\n', depUTC);
fprintf('Arrival UTC               : %s\n', arrUTC);
fprintf('Time of flight            : %.2f days\n', tofDays);
fprintf('Lambert position error    : %.6e km\n', posErr);
fprintf('Lambert velocity error    : %.6e km/s\n', velErr);
fprintf('Launch C3                 : %.6f km^2/s^2\n', C3);
fprintf('Departure v_inf           : %.6f km/s\n', vInfDep);
fprintf('Arrival   v_inf           : %.6f km/s\n', vInfArr);
fprintf('Earth departure dv        : %.6f km/s\n', dvDep);
fprintf('Mars  arrival dv          : %.6f km/s\n', dvArr);
fprintf('Total mission dv          : %.6f km/s\n', dvTot);
fprintf('Earth orbit umbra frac    : %.4f\n', mean(isUmbra));
fprintf('Earth access intervals    : %d\n', numel(intervals));
fprintf('============================================================\n');
fprintf('Done.\n');