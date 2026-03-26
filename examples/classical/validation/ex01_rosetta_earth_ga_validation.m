clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/astroToolbox')
startup

% ==========================================================
% Rosetta 1st Earth gravity-assist validation
% Explicit kernel loading, no local helper functions
% ==========================================================

% ----------------------------------------------------------
% Reset SPICE and load only the needed Rosetta kernels
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
% Time window around first Earth flyby
% ----------------------------------------------------------
tStartUTC = '2005-03-04 00:00:00';
tFlybyUTC = '2005-03-04 10:00:00';
tEndUTC   = '2005-03-11 00:00:00';

turnSign = +1;
altitude = 3000;   % km above Earth surface, first test value

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
    'VENUS BARYCENTER', ...
    'EARTH BARYCENTER', ...
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
fprintf('Rosetta 1st Earth gravity-assist validation\n');
fprintf('============================================================\n');
fprintf('Kernel root           : %s\n', kernelRoot);
fprintf('Spacecraft target     : %s\n', scTarget);
fprintf('Inbound start         : %s\n', tStartUTC);
fprintf('Flyby epoch           : %s\n', tFlybyUTC);
fprintf('Outbound end          : %s\n', tEndUTC);
fprintf('Chosen Earth altitude : %.3f km\n', altitude);
fprintf('Chosen turn sign      : %d\n', turnSign);

% ----------------------------------------------------------
% Initial state from SPICE
% ----------------------------------------------------------
state0 = cspice_spkezr(scTarget, etStart, 'J2000', 'NONE', 'SUN');
x0 = state0(:);

% ----------------------------------------------------------
% Inbound propagation
% ----------------------------------------------------------
fprintf('\nPropagating inbound arc...\n');
outIn = astro.propagators.propagateHighFidelity(x0, etStart, tofIn, configHF);

xMinus = outIn.x(end,:).';
rMinus = xMinus(1:3);
vMinus = xMinus(4:6);

stateEarth = cspice_spkezr('EARTH', etFlyby, 'J2000', 'NONE', 'SUN');
vEarth = stateEarth(4:6);

vInfIn = vMinus - vEarth;

fprintf('Incoming |v_inf^-|    : %.6f km/s\n', norm(vInfIn));

% ----------------------------------------------------------
% Flyby map using existing planar flybyTurn
% ----------------------------------------------------------
fb = astro.mga.flybyTurn(vInfIn, earth.mu, rpEarth, turnSign);
vPlus = vEarth + fb.vInfOut;

fprintf('Model turn angle      : %.6f deg\n', fb.deltaDeg);
fprintf('Outgoing |v_inf^+|    : %.6f km/s\n', norm(fb.vInfOut));

% ----------------------------------------------------------
% Outbound propagation
% ----------------------------------------------------------
fprintf('Propagating outbound arc...\n');
xFlybyPlus = [rMinus; vPlus];
outOut = astro.propagators.propagateHighFidelity(xFlybyPlus, etFlyby, tofOut, configHF);

% ----------------------------------------------------------
% Stitch reconstructed trajectory
% ----------------------------------------------------------
xRecon = [outIn.x; outOut.x(2:end,:)];
tRecon = [etStart + outIn.t; etFlyby + outOut.t(2:end)];

% ----------------------------------------------------------
% SPICE truth at same epochs
% ----------------------------------------------------------
fprintf('Sampling SPICE truth...\n');
xTruth = zeros(numel(tRecon), 6);

for k = 1:numel(tRecon)
    stateK = cspice_spkezr(scTarget, tRecon(k), 'J2000', 'NONE', 'SUN');
    xTruth(k,:) = stateK(:).';
end

% ----------------------------------------------------------
% Errors
% ----------------------------------------------------------
posErr = vecnorm(xRecon(:,1:3) - xTruth(:,1:3), 2, 2);
velErr = vecnorm(xRecon(:,4:6) - xTruth(:,4:6), 2, 2);

fprintf('\n============================================================\n');
fprintf('Reconstruction errors against SPICE\n');
fprintf('============================================================\n');
fprintf('Final position error : %.3f km\n', posErr(end));
fprintf('Final velocity error : %.6f km/s\n', velErr(end));
fprintf('Max position error   : %.3f km\n', max(posErr));
fprintf('Max velocity error   : %.6f km/s\n', max(velErr));
fprintf('RMS position error   : %.3f km\n', sqrt(mean(posErr.^2)));
fprintf('RMS velocity error   : %.6f km/s\n', sqrt(mean(velErr.^2)));

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
    'DisplayName', 'Rosetta SPICE truth');

plot(xRecon(:,1), xRecon(:,2), 'b-', ...
    'LineWidth', 1.8, ...
    'DisplayName', 'Reconstructed with Earth GA');

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
title('Rosetta 1st Earth Gravity Assist: Reconstruction vs SPICE', ...
    'FontSize', 15, 'FontWeight', 'bold');

legend('Location', 'best');

fprintf('\nDone.\n');