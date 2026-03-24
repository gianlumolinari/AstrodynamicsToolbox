clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/astroToolbox')
startup

addpath(genpath('/Users/gianlucamolinari/Desktop/astroToolbox/external/mice'))
rehash

fprintf('\n');
fprintf('============================================================\n');
fprintf('HERA SPICE Trajectory Plot\n');
fprintf('============================================================\n');

kernelRoot = '/Users/gianlucamolinari/Desktop/astroToolbox/data/spice/HERA/kernels';
mkDir = fullfile(kernelRoot, 'mk');

mkFile = fullfile(mkDir, 'hera_crema_2_1.tm');
if exist(mkFile, 'file') ~= 2
    error('Meta-kernel not found: %s', mkFile);
end

cspice_kclear

oldDir = pwd;
cleanupObj = onCleanup(@() cd(oldDir)); 

cd(mkDir)
cspice_furnsh(mkFile)

heraID = cspice_bodn2c('HERA');
fprintf('HERA NAIF ID: %d\n', heraID);

% Only use SPKs that actually contain HERA = -91
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

    fprintf('Using SPK: %s\n', spkNames{i});
end

if isempty(coverAll)
    cspice_kclear
    error('No HERA coverage found in selected SPKs.');
end

nIntervals = cspice_wncard(coverAll);
fprintf('Coverage intervals found: %d\n', nIntervals);

[leftET, rightET] = cspice_wnfetd(coverAll, 1);

fprintf('Coverage start : %s\n', cspice_et2utc(leftET,  'C', 3));
fprintf('Coverage stop  : %s\n', cspice_et2utc(rightET, 'C', 3));
fprintf('Coverage days  : %.6f\n', (rightET-leftET)/86400);

nSamples = 1500;
etGridAll = linspace(leftET, rightET, nSamples);

rHera  = [];
rEarth = [];
etGrid = [];

for k = 1:nSamples
    et = etGridAll(k);

    try
        stH = cspice_spkezr('HERA',  et, 'J2000', 'NONE', 'SUN');
        stE = cspice_spkezr('EARTH', et, 'J2000', 'NONE', 'SUN');

        rHera(end+1,:)  = stH(1:3).'; %#ok<AGROW>
        rEarth(end+1,:) = stE(1:3).'; %#ok<AGROW>
        etGrid(end+1,1) = et; %#ok<AGROW>
    catch
        % Skip epochs where the full HERA->SUN state chain is unavailable
    end
end

fprintf('Valid sampled epochs for HERA wrt SUN: %d / %d\n', numel(etGrid), nSamples);

if isempty(etGrid)
    cspice_kclear
    error('No valid HERA heliocentric states were found over the sampled interval.');
end

fprintf('First valid epoch: %s\n', cspice_et2utc(etGrid(1),   'C', 3));
fprintf('Last valid epoch : %s\n', cspice_et2utc(etGrid(end), 'C', 3));

figure('Color','w');
hold on
grid on
box on
axis equal

plot(rEarth(:,1), rEarth(:,2), 'b-', 'LineWidth', 1.2);
plot(rHera(:,1),  rHera(:,2),  'r-', 'LineWidth', 1.8);

plot(rEarth(1,1), rEarth(1,2), 'bo', 'MarkerSize', 7, 'MarkerFaceColor', 'b');
plot(rEarth(end,1), rEarth(end,2), 'bs', 'MarkerSize', 7, 'LineWidth', 1.2);

plot(rHera(1,1), rHera(1,2), 'ro', 'MarkerSize', 7, 'MarkerFaceColor', 'r');
plot(rHera(end,1), rHera(end,2), 'rs', 'MarkerSize', 7, 'LineWidth', 1.2);

plot(0, 0, 'yo', 'MarkerSize', 10, 'LineWidth', 1.5, 'MarkerFaceColor', 'y');

xlabel('x [km]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('y [km]', 'FontSize', 13, 'FontWeight', 'bold');
title('HERA SPICE Truth Trajectory - Top View (Heliocentric J2000)', ...
    'FontSize', 15, 'FontWeight', 'bold');

legend('Earth', 'HERA', 'Earth start', 'Earth end', ...
       'HERA start', 'HERA end', 'Sun', 'Location', 'best');

view(2)

fprintf('\nDone.\n');

cspice_kclear