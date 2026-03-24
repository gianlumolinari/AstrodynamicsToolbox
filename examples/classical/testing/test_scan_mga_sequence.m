clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/AstroToolbox')
startup
astro.ephem.loadSpiceKernels(fullfile(pwd, 'data', 'spice'));

sun = astro.bodies.getBody('sun');

% ----------------------------------------------------------
% Choose sequence
% ----------------------------------------------------------
sequence = {'earth','earth','jupiter'};
% sequence = {'earth','venus','earth','earth','jupiter'};

% ----------------------------------------------------------
% Departure window
% ----------------------------------------------------------
depDates = datetime(2025,1,1,0,0,0):days(30):datetime(2028,1,1,0,0,0);
depDatesUTC = cell(numel(depDates),1);
for k = 1:numel(depDates)
    depDatesUTC{k} = datestr(depDates(k), 'yyyy-mm-dd HH:MM:SS');
end

% ----------------------------------------------------------
% TOF grids per leg [days]
% Adjust these depending on sequence length
% ----------------------------------------------------------
if numel(sequence) == 3
    % Earth -> Earth -> Jupiter
    tofGridDays = {
        200:40:600
        500:80:1800
    };

    bodyMap.earth = astro.bodies.getBody('earth');
    minAltMap.earth = 300;

elseif numel(sequence) == 5
    % Earth -> Venus -> Earth -> Earth -> Jupiter
    tofGridDays = {
        80:40:300
        80:40:300
        120:40:500
        500:100:1800
    };

    bodyMap.venus = astro.bodies.getBody('venus');
    bodyMap.earth = astro.bodies.getBody('earth');

    minAltMap.venus = 300;
    minAltMap.earth = 300;
else
    error('This example currently supports 3-body or 5-body sequences.');
end

fprintf('Running coarse MGA scan...\n');
results = astro.mga.scanMGAWindow(sequence, depDatesUTC, tofGridDays, sun.mu, bodyMap, minAltMap, 'spice');

if isempty(results)
    error('No valid scan points were generated.');
end

nShow = min(10, numel(results));

fprintf('\nTop %d candidates by objective:\n', nShow);
fprintf('------------------------------------------\n');

for k = 1:nShow
    r = results(k);
    fprintf('\nCandidate %d\n', k);
    fprintf('  Departure UTC           : %s\n', r.depUTC);
    fprintf('  TOFs [days]             : ');
    fprintf('%.2f ', r.tofDays);
    fprintf('\n');
    fprintf('  Objective               : %.6f\n', r.objective);
    fprintf('  Max |v_inf| mismatch    : %.6f km/s\n', r.maxMagMismatch);
    fprintf('  Max turn violation      : %.6f deg\n', r.maxTurnViolationDeg);
    fprintf('  All flybys feasible     : %d\n', r.allFlybysFeasible);
end

% ----------------------------------------------------------
% Plot objective of top candidates
% ----------------------------------------------------------
objVals = [results(1:nShow).objective];

figure('Color','w');
plot(1:nShow, objVals, 'o-', 'LineWidth', 1.3);
grid on
xlabel('Rank', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Objective', 'FontSize', 13, 'FontWeight', 'bold');
title('Top MGA Scan Candidates', 'FontSize', 15, 'FontWeight', 'bold');