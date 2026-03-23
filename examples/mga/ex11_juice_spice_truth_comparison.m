clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/AstroToolbox')
startup

% Load generic planetary kernels first
astro.ephem.loadSpiceKernels(fullfile(pwd, 'data', 'spice'));

% ==========================================================
% JUICE SPICE truth comparison
%
% This script:
%   1) loads JUICE mission kernels from data/spice/juice
%   2) reconstructs a reference Lambert-leg sequence
%   3) reads the real JUICE spacecraft trajectory from SPICE
%   4) overlays reconstructed and real trajectories
%   5) compares encounter states at selected epochs
% ==========================================================

sun = astro.bodies.getBody('sun');

% ----------------------------------------------------------
% Reference heliocentric sequence used for reconstruction
% ----------------------------------------------------------
sequence = {'earth','earth','venus','earth','earth','jupiter'};

dates = {
    '2023-04-14 12:14:00'
    '2024-08-20 21:56:00'
    '2025-08-31 00:00:00'
    '2026-09-01 00:00:00'
    '2029-01-01 00:00:00'
    '2031-07-01 00:00:00'
};

planetNames = {'Earth','Earth (LEGA)','Venus','Earth','Earth','Jupiter'};

fprintf('\n');
fprintf('============================================================\n');
fprintf('JUICE SPICE Truth Comparison\n');
fprintf('============================================================\n');

% ----------------------------------------------------------
% Load JUICE mission kernels
% ----------------------------------------------------------
juiceKernelDir = fullfile(pwd, 'data', 'spice', 'juice');

if exist(juiceKernelDir, 'dir') ~= 7
    error(['JUICE kernel directory not found: ', juiceKernelDir]);
end

fprintf('\nLoading JUICE mission kernels from:\n  %s\n', juiceKernelDir);
localLoadJuiceKernels(juiceKernelDir);

% ----------------------------------------------------------
% Try to resolve JUICE spacecraft target name
% ----------------------------------------------------------
candidateTargets = {'JUICE', 'JUICE SPACECRAFT', '-28'};
juiceTarget = '';

testEpoch = dates{1};
for k = 1:numel(candidateTargets)
    try
        tmp = astro.ephem.getSpiceState(candidateTargets{k}, testEpoch, 'SUN', 'J2000', 'NONE');
        juiceTarget = tmp.target;
        break
    catch
    end
end

if isempty(juiceTarget)
    fprintf('\nCould not resolve JUICE spacecraft target automatically.\n');
    fprintf('Available loaded kernels may use a different spacecraft ID/name.\n');
    fprintf('You may need to inspect the mission FK/SPK naming.\n');
    error('JUICE target name resolution failed.');
else
    fprintf('\nResolved JUICE spacecraft target as: %s\n', juiceTarget);
end

% ----------------------------------------------------------
% Reconstruct reference Lambert-leg sequence
% ----------------------------------------------------------
traj = astro.mga.propagateMGA(sequence, dates, sun.mu, 'spice');

% ----------------------------------------------------------
% Sample real JUICE spacecraft trajectory from SPICE
% ----------------------------------------------------------
t0 = datetime(dates{1}, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
tf = datetime(dates{end}, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');

nSamples = 1200;
tHist = linspace(posixtime(t0), posixtime(tf), nSamples);
tHist = datetime(tHist, 'ConvertFrom', 'posixtime');

rJuice = zeros(nSamples, 3);

fprintf('\nSampling real JUICE trajectory from SPICE...\n');
for k = 1:nSamples
    epochUTC = datestr(tHist(k), 'yyyy-mm-dd HH:MM:SS');
    sc = astro.ephem.getSpiceState(juiceTarget, epochUTC, 'SUN', 'J2000', 'NONE');
    rJuice(k,:) = sc.r(:).';
end

% ----------------------------------------------------------
% Build reconstructed Lambert trajectory samples
% ----------------------------------------------------------
rRecon = [];
for k = 1:numel(traj.legs)
    x0Leg = [traj.legs(k).r1; traj.legs(k).v1];
    tspan = [0, traj.legs(k).tofSec];

    optsProp.RelTol = 1e-11;
    optsProp.AbsTol = 1e-11;
    optsProp.Solver = 'ode113';

    out = astro.propagators.propagate( ...
        @(t,x) astro.propagators.eomTwoBody(t, x, sun.mu), ...
        tspan, x0Leg, optsProp);

    rRecon = [rRecon; out.x(:,1:3)]; %#ok<AGROW>
end

% ----------------------------------------------------------
% Plot comparison
% ----------------------------------------------------------
figure('Color','w');
hold on
axis equal
grid on

hSun   = plot(0, 0, 'ko', 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', 'Sun');
hRecon = plot(rRecon(:,1), rRecon(:,2), 'b-', 'LineWidth', 1.3, 'DisplayName', 'Lambert reconstruction');
hTruth = plot(rJuice(:,1), rJuice(:,2), 'g--', 'LineWidth', 1.4, 'DisplayName', 'JUICE SPICE truth');

% Encounter markers
hTruthPts = plot(nan, nan, 'go', 'MarkerSize', 7, 'LineWidth', 1.2, 'DisplayName', 'JUICE encounter states');
hReconPts = plot(nan, nan, 'bs', 'MarkerSize', 7, 'LineWidth', 1.2, 'DisplayName', 'Reconstructed encounter states');

for k = 1:numel(dates)
    epochUTC = dates{k};
    sc = astro.ephem.getSpiceState(juiceTarget, epochUTC, 'SUN', 'J2000', 'NONE');
    plot(sc.r(1), sc.r(2), 'go', 'MarkerSize', 7, 'LineWidth', 1.2);

    if k == 1
        rRef = traj.legs(1).r1;
    elseif k == numel(dates)
        rRef = traj.legs(end).r2;
    else
        rRef = traj.legs(k-1).r2;
    end

    plot(rRef(1), rRef(2), 'bs', 'MarkerSize', 7, 'LineWidth', 1.2);

    % Planet/encounter labels near SPICE-truth states
    text(sc.r(1), sc.r(2), ['  ' planetNames{k}], ...
        'FontSize', 10, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
end

xlabel('x [km]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('y [km]', 'FontSize', 13, 'FontWeight', 'bold');
title('JUICE: SPICE Truth vs Lambert Reconstruction', ...
    'FontSize', 15, 'FontWeight', 'bold');
legend([hSun, hRecon, hTruth, hTruthPts, hReconPts], ...
       'Location', 'best');

% ----------------------------------------------------------
% Epoch-by-epoch position comparison
% ----------------------------------------------------------
fprintf('\nEncounter state comparison:\n');
fprintf('------------------------------------------------------------\n');
for k = 1:numel(dates)
    epochUTC = dates{k};
    sc = astro.ephem.getSpiceState(juiceTarget, epochUTC, 'SUN', 'J2000', 'NONE');

    if k == 1
        rRef = traj.legs(1).r1;
        vRef = traj.legs(1).v1;
    elseif k == numel(dates)
        rRef = traj.legs(end).r2;
        vRef = traj.legs(end).v2;
    else
        rRef = traj.legs(k-1).r2;
        vRef = traj.legs(k-1).v2;
    end

    posErr = norm(sc.r - rRef);
    velErr = norm(sc.v - vRef);

    fprintf('%-19s  pos err = %12.3f km   vel err = %10.6f km/s\n', ...
        epochUTC, posErr, velErr);
end

fprintf('\nDone.\n');

function localLoadJuiceKernels(kernelDir)
    files = dir(fullfile(kernelDir, '**', '*'));
    for k = 1:numel(files)
        if files(k).isdir
            continue
        end
        [~,~,ext] = fileparts(files(k).name);
        if any(strcmpi(ext, {'.bsp','.bc','.tf','.tls','.tpc','.tsc','.tm','.bpc'}))
            try
                cspice_furnsh(fullfile(files(k).folder, files(k).name));
            catch
            end
        end
    end
end