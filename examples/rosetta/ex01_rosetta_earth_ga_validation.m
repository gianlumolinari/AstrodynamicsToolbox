clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/astroToolbox')
startup

% ==========================================================
% Rosetta Earth flyby:
% SPICE trajectory vs manually reconstructed MGA trajectory
% ==========================================================

% ----------------------------------------------------------
% Reset SPICE and load Rosetta kernels
% ----------------------------------------------------------
cspice_kclear;

kernelRoot = fullfile(pwd, 'data', 'spice', 'ROSETTA', 'kernels');

cspice_furnsh(fullfile(kernelRoot, 'lsk', 'NAIF0011.TLS'));
cspice_furnsh(fullfile(kernelRoot, 'pck', 'PCK00010.TPC'));
cspice_furnsh(fullfile(kernelRoot, 'spk', 'DE405.BSP'));
cspice_furnsh(fullfile(kernelRoot, 'spk', 'ORHR___________T19_00122.BSP'));

% ----------------------------------------------------------
% Spacecraft target
% ----------------------------------------------------------
scTarget = '-226';

% ----------------------------------------------------------
% Manual flyby setup
% ----------------------------------------------------------
tStartUTC = '2005-03-04 10:00:00';
tFlybyUTC = '2005-03-04 22:09:06';
tEndUTC   = '2005-03-06 10:00:00';

% Choose ONE of the two models below:
use3DGA = true;     % true -> gravityAssist3D, false -> flybyTurn

% Manual flyby parameters
altitude = 3000;    % km above Earth surface
turnSign = -1;      % used only if use3DGA = false
thetaB   = deg2rad(-110);   % used only if use3DGA = true

% ----------------------------------------------------------
% Bodies
% ----------------------------------------------------------
mercury = astro.bodies.getBody('mercury');
venus   = astro.bodies.getBody('venus');
earth   = astro.bodies.getBody('earth');
moon    = astro.bodies.getBody('moon');
mars    = astro.bodies.getBody('mars');
jupiter = astro.bodies.getBody('jupiter');
saturn  = astro.bodies.getBody('saturn');

rpEarth = earth.radius + altitude;

% ----------------------------------------------------------
% HF propagator config
% ----------------------------------------------------------
configHF = struct();
configHF.RelTol   = 1e-12;
configHF.AbsTol   = 1e-12;
configHF.Cr       = 1.2;
configHF.A_over_m = 0.0;

configHF.bodyNames = { ...
    'MERCURY BARYCENTER', ...
    'VENUS', ...
    'EARTH', ...
    'MOON', ...
    'MARS BARYCENTER', ...
    'JUPITER BARYCENTER', ...
    'SATURN BARYCENTER'};

configHF.muBodies = [ ...
    mercury.mu, ...
    venus.mu, ...
    earth.mu, ...
    moon.mu, ...
    mars.mu, ...
    jupiter.mu, ...
    saturn.mu];

% ----------------------------------------------------------
% Epochs
% ----------------------------------------------------------
etStart = cspice_str2et(tStartUTC);
etFlyby = cspice_str2et(tFlybyUTC);
etEnd   = cspice_str2et(tEndUTC);

tofIn  = etFlyby - etStart;
tofOut = etEnd   - etFlyby;

if tofIn <= 0 || tofOut <= 0
    error('Check tStartUTC, tFlybyUTC, tEndUTC chronology.');
end

fprintf('\n============================================================\n');
fprintf('Rosetta Earth flyby: SPICE vs MGA reconstruction\n');
fprintf('============================================================\n');
fprintf('Spacecraft target     : %s\n', scTarget);
fprintf('Start epoch           : %s\n', tStartUTC);
fprintf('Flyby epoch           : %s\n', tFlybyUTC);
fprintf('End epoch             : %s\n', tEndUTC);
fprintf('Altitude              : %.3f km\n', altitude);

if use3DGA
    fprintf('Flyby model           : gravityAssist3D\n');
    fprintf('thetaB                : %.3f deg\n', rad2deg(thetaB));
else
    fprintf('Flyby model           : flybyTurn\n');
    fprintf('turnSign              : %d\n', turnSign);
end

% ----------------------------------------------------------
% Initial state from SPICE
% ----------------------------------------------------------
state0 = cspice_spkezr(scTarget, etStart, 'J2000', 'NONE', 'SUN');
x0 = state0(:);

% ----------------------------------------------------------
% Inbound propagation
% ----------------------------------------------------------
outIn = astro.propagators.propagateHighFidelity(x0, etStart, tofIn, configHF);

xMinus = outIn.x(end,:).';
rMinus = xMinus(1:3);
vMinus = xMinus(4:6);

stateEarth = cspice_spkezr('EARTH', etFlyby, 'J2000', 'NONE', 'SUN');
vEarth = stateEarth(4:6);

vInfIn = vMinus - vEarth;

% ----------------------------------------------------------
% Apply manual gravity assist
% ----------------------------------------------------------
if use3DGA
    ga = astro.mga.gravityAssist3D('earth', rpEarth, thetaB, vInfIn, vEarth, etFlyby);
    vPlus = ga.vOutHelio;
    fprintf('Model turn angle      : %.6f deg\n', ga.deltaDeg);
    fprintf('|v_inf^-|             : %.6f km/s\n', ga.vInfMagIn);
    fprintf('|v_inf^+|             : %.6f km/s\n', ga.vInfMagOut);
else
    fb = astro.mga.flybyTurn(vInfIn, earth.mu, rpEarth, turnSign);
    vPlus = vEarth + fb.vInfOut;
    fprintf('Model turn angle      : %.6f deg\n', fb.deltaDeg);
    fprintf('|v_inf^-|             : %.6f km/s\n', fb.vInfMag);
    fprintf('|v_inf^+|             : %.6f km/s\n', norm(fb.vInfOut));
end

% ----------------------------------------------------------
% Outbound propagation
% ----------------------------------------------------------
xFlybyPlus = [rMinus; vPlus];
outOut = astro.propagators.propagateHighFidelity(xFlybyPlus, etFlyby, tofOut, configHF);

% ----------------------------------------------------------
% Stitch reconstructed trajectory
% ----------------------------------------------------------
xRecon = [outIn.x; outOut.x(2:end,:)];
tRecon = [etStart + outIn.t; etFlyby + outOut.t(2:end)];

% ----------------------------------------------------------
% Sample SPICE truth on same epochs
% ----------------------------------------------------------
xTruth = zeros(numel(tRecon), 6);

for k = 1:numel(tRecon)
    stateK = cspice_spkezr(scTarget, tRecon(k), 'J2000', 'NONE', 'SUN');
    xTruth(k,:) = stateK(:).';
end

% ----------------------------------------------------------
% Plot
% ----------------------------------------------------------
figure('Color','w');
hold on;
grid on;
axis equal;

plot(0, 0, 'ko', ...
    'MarkerSize', 10, ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Sun');

plot(xTruth(:,1), xTruth(:,2), 'k--', ...
    'LineWidth', 1.6, ...
    'DisplayName', 'Rosetta SPICE');

plot(xRecon(:,1), xRecon(:,2), 'b-', ...
    'LineWidth', 1.8, ...
    'DisplayName', 'MGA propagated');

plot(xTruth(1,1), xTruth(1,2), 'go', ...
    'MarkerSize', 7, ...
    'LineWidth', 1.2, ...
    'DisplayName', 'Start');

plot(xTruth(end,1), xTruth(end,2), 'rs', ...
    'MarkerSize', 7, ...
    'LineWidth', 1.2, ...
    'DisplayName', 'End');

stateFlybyTruth = cspice_spkezr(scTarget, etFlyby, 'J2000', 'NONE', 'SUN');
plot(stateFlybyTruth(1), stateFlybyTruth(2), 'mo', ...
    'MarkerSize', 7, ...
    'LineWidth', 1.2, ...
    'DisplayName', 'Flyby epoch');

xlabel('x [km]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('y [km]', 'FontSize', 13, 'FontWeight', 'bold');
title('Rosetta Earth Flyby: SPICE vs MGA Propagated Trajectory', ...
    'FontSize', 15, 'FontWeight', 'bold');

legend('Location', 'best');

fprintf('\nDone.\n');