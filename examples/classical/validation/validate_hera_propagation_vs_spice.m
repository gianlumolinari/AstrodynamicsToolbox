clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/astroToolbox')
startup

addpath(genpath('/Users/gianlucamolinari/Desktop/astroToolbox/external/mice'))
rehash

fprintf('\n');
fprintf('============================================================\n');
fprintf('HERA Piecewise High-Fidelity Propagation vs SPICE Validation\n');
fprintf('============================================================\n');

kernelRoot = '/Users/gianlucamolinari/Desktop/astroToolbox/data/spice/HERA/kernels';
mkDir      = fullfile(kernelRoot, 'mk');
mkFile     = fullfile(mkDir, 'hera_crema_2_1.tm');

if exist(mkFile, 'file') ~= 2
    error('Meta-kernel not found: %s', mkFile);
end

cspice_kclear

oldDir = pwd;
cleanupObj = onCleanup(@() cd(oldDir)); %#ok<NASGU>

cd(mkDir)
cspice_furnsh(mkFile)

% ----------------------------------------------------------
% Resolve HERA and build candidate heliocentric interval
% ----------------------------------------------------------
heraID = cspice_bodn2c('HERA');
fprintf('HERA NAIF ID: %d\n', heraID);

spkNames = { ...
    'hera_fcp_000067_241007_261104_v01.bsp', ...
    'hera_flp_000005_241007_270430_v01.bsp', ...
    'hera_sc_LPO_EMA_2024c_v01.bsp', ...
    'hera_sc_PO_LPO_EMA_2024_v05.bsp', ...
    'hera_sc_crema_2_1_ECP_PDP_DCP_261125_270303_v01.bsp', ...
    'hera_sc_crema_2_1_LPO_241007_261202_v01.bsp'};

room = 10000;
coverAll = [];

for i = 1:numel(spkNames)
    spkPath = fullfile(kernelRoot, 'spk', spkNames{i});

    if exist(spkPath, 'file') ~= 2
        fprintf('Missing SPK: %s\n', spkNames{i});
        continue
    end

    coverSeg = cspice_spkcov(spkPath, heraID, room);

    if isempty(coverAll)
        coverAll = coverSeg;
    else
        coverAll = cspice_wnunid(coverAll, coverSeg);
    end
end

if isempty(coverAll)
    cspice_kclear
    error('No HERA object coverage found.');
end

[leftET, rightET] = cspice_wnfetd(coverAll, 1);

fprintf('Candidate coverage start : %s\n', cspice_et2utc(leftET,  'C', 3));
fprintf('Candidate coverage stop  : %s\n', cspice_et2utc(rightET, 'C', 3));

% ----------------------------------------------------------
% Build a valid SPICE truth interval for HERA wrt SUN
% ----------------------------------------------------------
nProbe  = 4000;
etProbe = linspace(leftET, rightET, nProbe);

validMask = false(size(etProbe));

for k = 1:nProbe
    try
        cspice_spkezr('HERA', etProbe(k), 'J2000', 'NONE', 'SUN');
        validMask(k) = true;
    catch
        validMask(k) = false;
    end
end

if ~any(validMask)
    cspice_kclear
    error('No valid HERA heliocentric states found.');
end

idxFirst = find(validMask, 1, 'first');
idxLast  = find(validMask, 1, 'last');

leftValidET  = etProbe(idxFirst);
rightValidET = etProbe(idxLast);

fprintf('Valid HERA->SUN start    : %s\n', cspice_et2utc(leftValidET,  'C', 3));
fprintf('Valid HERA->SUN stop     : %s\n', cspice_et2utc(rightValidET, 'C', 3));
fprintf('Valid interval span [d]  : %.6f\n', (rightValidET-leftValidET)/86400);

% ----------------------------------------------------------
% Sample SPICE truth over full valid interval
% ----------------------------------------------------------
nSamples = 2000;
etGrid   = linspace(leftValidET, rightValidET, nSamples);

Xtruth = zeros(nSamples, 6);
rEarth = zeros(nSamples, 3);

for k = 1:nSamples
    stH = cspice_spkezr('HERA',  etGrid(k), 'J2000', 'NONE', 'SUN');
    stE = cspice_spkezr('EARTH', etGrid(k), 'J2000', 'NONE', 'SUN');

    Xtruth(k,:) = stH(:).';
    rEarth(k,:) = stE(1:3).';
end

x0 = Xtruth(1,:).';

fprintf('\nInitial SPICE state used for propagation:\n');
disp(x0)

% ----------------------------------------------------------
% High-fidelity propagator configuration
% ----------------------------------------------------------
sun     = astro.bodies.getBody('sun');
earth   = astro.bodies.getBody('earth');
moon    = astro.bodies.getBody('moon');
venus   = astro.bodies.getBody('venus');
mars    = astro.bodies.getBody('mars');
jupiter = astro.bodies.getBody('jupiter');
saturn  = astro.bodies.getBody('saturn');

config = struct();
config.muSun = sun.mu;

config.bodyNames = { ...
    'VENUS BARYCENTER', ...
    'EARTH', ...
    'MOON', ...
    'MARS BARYCENTER', ...
    'JUPITER BARYCENTER', ...
    'SATURN BARYCENTER'};

config.muBodies = [ ...
    venus.mu, ...
    earth.mu, ...
    moon.mu, ...
    mars.mu, ...
    jupiter.mu, ...
    saturn.mu];

config.useSRP   = true;
config.Cr       = 1.3;
config.A_over_m = 0.01;

config.RelTol = 1e-12;
config.AbsTol = 1e-12;

% ----------------------------------------------------------
% Piecewise propagation settings
% ----------------------------------------------------------
segmentDays = 20;                  % chunk size
segmentSec  = segmentDays * 86400; % [s]

tTotal = etGrid - leftValidET;
tfinal = tTotal(end);

% Segment boundaries in time-since-start coordinates
tBreaks = 0:segmentSec:tfinal;
if tBreaks(end) < tfinal
    tBreaks(end+1) = tfinal;
end

fprintf('\nPiecewise propagation:\n');
fprintf('  Segment length [days] : %.3f\n', segmentDays);
fprintf('  Number of segments    : %d\n', numel(tBreaks)-1);

Xmodel = zeros(nSamples, 6);

xInitSeg = x0;

for iSeg = 1:numel(tBreaks)-1
    tStartSeg = tBreaks(iSeg);
    tEndSeg   = tBreaks(iSeg+1);

    fprintf('  Segment %3d: %10.3f -> %10.3f days\n', ...
        iSeg, tStartSeg/86400, tEndSeg/86400);

    % Epoch at start of current segment
    etStartSeg = leftValidET + tStartSeg;

    % Update config epoch for this segment
    config.et0 = etStartSeg;

    % Propagate current segment
    outSeg = astro.propagators.propagateHighFidelity( ...
        xInitSeg, etStartSeg, tEndSeg - tStartSeg, config);

    % Global sample indices belonging to this segment
    if iSeg < numel(tBreaks)-1
        idxSeg = find(tTotal >= tStartSeg & tTotal < tEndSeg);
    else
        idxSeg = find(tTotal >= tStartSeg & tTotal <= tEndSeg);
    end

    % Interpolate segment solution onto global SPICE epochs
    tLocal = tTotal(idxSeg) - tStartSeg;

    for j = 1:6
        Xmodel(idxSeg, j) = interp1(outSeg.t, outSeg.x(:,j), tLocal, 'pchip');
    end

    % Carry final propagated state into next segment
    xInitSeg = outSeg.x(end,:).';
end

% ----------------------------------------------------------
% Compute errors
% ----------------------------------------------------------
dr = Xmodel(:,1:3) - Xtruth(:,1:3);
dv = Xmodel(:,4:6) - Xtruth(:,4:6);

rErr = vecnorm(dr, 2, 2);
vErr = vecnorm(dv, 2, 2);

fprintf('\nError summary against SPICE truth:\n');
fprintf('  Initial position error [km]   : %.6e\n', rErr(1));
fprintf('  Final position error [km]     : %.6e\n', rErr(end));
fprintf('  Max position error [km]       : %.6e\n', max(rErr));
fprintf('  RMS position error [km]       : %.6e\n', sqrt(mean(rErr.^2)));

fprintf('  Initial velocity error [km/s] : %.6e\n', vErr(1));
fprintf('  Final velocity error [km/s]   : %.6e\n', vErr(end));
fprintf('  Max velocity error [km/s]     : %.6e\n', max(vErr));
fprintf('  RMS velocity error [km/s]     : %.6e\n', sqrt(mean(vErr.^2)));

[~, idxRmax] = max(rErr);
[~, idxVmax] = max(vErr);

fprintf('\nAdditional diagnostics:\n');
fprintf('  Max position error epoch [UTC] : %s\n', cspice_et2utc(etGrid(idxRmax), 'C', 3));
fprintf('  Max velocity error epoch [UTC] : %s\n', cspice_et2utc(etGrid(idxVmax), 'C', 3));

% ----------------------------------------------------------
% Graphical comparison
% ----------------------------------------------------------
figure('Color', 'w');
hold on
grid on
box on
axis equal

plot(rEarth(:,1), rEarth(:,2), 'b-',  'LineWidth', 1.2);
plot(Xtruth(:,1), Xtruth(:,2), 'k-',  'LineWidth', 1.8);
plot(Xmodel(:,1), Xmodel(:,2), 'r--', 'LineWidth', 1.6);

plot(0, 0, 'yo', 'MarkerSize', 10, 'LineWidth', 1.4, 'MarkerFaceColor', 'y');

plot(Xtruth(1,1),   Xtruth(1,2),   'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'k');
plot(Xtruth(end,1), Xtruth(end,2), 'ks', 'MarkerSize', 7, 'LineWidth', 1.2);

xlabel('x [km]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('y [km]', 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('HERA: SPICE vs Piecewise High-Fidelity Propagation (%g-day segments)', segmentDays), ...
      'FontSize', 15, 'FontWeight', 'bold');
legend('Earth (SPICE)', 'HERA truth (SPICE)', 'HERA propagated (toolbox)', ...
       'Sun', 'Truth start', 'Truth end', 'Location', 'best');
view(2)

figure('Color', 'w');
hold on
grid on
box on
plot(tTotal/86400, rErr, 'k-', 'LineWidth', 1.5);
xlabel('Time since initial epoch [days]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Position error [km]', 'FontSize', 13, 'FontWeight', 'bold');
title('Position Error vs SPICE Truth', 'FontSize', 15, 'FontWeight', 'bold');

figure('Color', 'w');
hold on
grid on
box on
plot(tTotal/86400, vErr, 'k-', 'LineWidth', 1.5);
xlabel('Time since initial epoch [days]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Velocity error [km/s]', 'FontSize', 13, 'FontWeight', 'bold');
title('Velocity Error vs SPICE Truth', 'FontSize', 15, 'FontWeight', 'bold');

fprintf('\nValidation completed.\n');

cspice_kclear