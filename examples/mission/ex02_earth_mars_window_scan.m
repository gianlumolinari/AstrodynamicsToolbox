clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/astroToolbox')
startup

% ==========================================================
% Earth-Mars porkchop-style window scan using Horizons + Izzo
% ==========================================================

sun = astro.bodies.getBody('sun');

% ----------------------------------------------------------
% User settings
% ----------------------------------------------------------
depStart = datetime('2026-09-01 00:00:00', ...
    'InputFormat', 'yyyy-MM-dd HH:mm:ss', 'TimeZone', 'UTC');
depEnd   = datetime('2026-12-01 00:00:00', ...
    'InputFormat', 'yyyy-MM-dd HH:mm:ss', 'TimeZone', 'UTC');

arrStart = datetime('2027-05-01 00:00:00', ...
    'InputFormat', 'yyyy-MM-dd HH:mm:ss', 'TimeZone', 'UTC');
arrEnd   = datetime('2027-11-01 00:00:00', ...
    'InputFormat', 'yyyy-MM-dd HH:mm:ss', 'TimeZone', 'UTC');

% Finer resolution for smoother porkchops
depStepDays = 2;
arrStepDays = 2;

minTOFdays = 60;

depDates = depStart:days(depStepDays):depEnd;
arrDates = arrStart:days(arrStepDays):arrEnd;

nDep = numel(depDates);
nArr = numel(arrDates);

% ----------------------------------------------------------
% Preallocate result grids
% Rows -> arrival dates
% Cols -> departure dates
% ----------------------------------------------------------
C3Grid      = NaN(nArr, nDep);
vInfDepGrid = NaN(nArr, nDep);
vInfArrGrid = NaN(nArr, nDep);
TOFGrid     = NaN(nArr, nDep);

% ----------------------------------------------------------
% Pre-fetch ephemerides
% ----------------------------------------------------------
earthStates = cell(1, nDep);
for j = 1:nDep
    epochStr = datestr(depDates(j), 'yyyy-mm-dd HH:MM:SS');
    earthStates{j} = astro.ephem.getHorizonsState('399', epochStr, '500@10');
end

marsStates = cell(1, nArr);
for i = 1:nArr
    epochStr = datestr(arrDates(i), 'yyyy-mm-dd HH:MM:SS');
    marsStates{i} = astro.ephem.getHorizonsState('499', epochStr, '500@10');
end

% ----------------------------------------------------------
% Main scan loop
% ----------------------------------------------------------
for i = 1:nArr
    for j = 1:nDep
        
        tof = seconds(arrDates(i) - depDates(j));
        tofDays = tof / 86400;
        
        if tofDays <= minTOFdays
            continue
        end
        
        r1 = earthStates{j}.r;
        vP1 = earthStates{j}.v;
        
        r2 = marsStates{i}.r;
        vP2 = marsStates{i}.v;
        
        try
            sol = astro.lambert.solveIzzo(r1, r2, tof, sun.mu, false);
            
            if ~sol.converged
                continue
            end
            
            vInfDepVec = astro.maneuvers.hyperbolicExcess(sol.v1, vP1);
            vInfArrVec = astro.maneuvers.hyperbolicExcess(sol.v2, vP2);
            
            vInfDep = norm(vInfDepVec);
            vInfArr = norm(vInfArrVec);
            c3 = astro.maneuvers.C3(vInfDep);
            
            C3Grid(i,j)      = c3;
            vInfDepGrid(i,j) = vInfDep;
            vInfArrGrid(i,j) = vInfArr;
            TOFGrid(i,j)     = tofDays;
            
        catch
            continue
        end
    end
end

% ----------------------------------------------------------
% Best points
% ----------------------------------------------------------
validC3 = C3Grid;
validDep = vInfDepGrid;
validArr = vInfArrGrid;

[minC3, idxC3] = min(validC3(:), [], 'omitnan');
[minDepVinf, idxDep] = min(validDep(:), [], 'omitnan');
[minArrVinf, idxArr] = min(validArr(:), [], 'omitnan');

if ~isnan(minC3)
    [iC3, jC3] = ind2sub(size(C3Grid), idxC3);
    fprintf('Minimum C3 found = %.6f km^2/s^2\n', minC3);
    fprintf('  Departure: %s\n', datestr(depDates(jC3)));
    fprintf('  Arrival  : %s\n', datestr(arrDates(iC3)));
    fprintf('  TOF      : %.2f days\n\n', TOFGrid(iC3,jC3));
else
    iC3 = [];
    jC3 = [];
end

if ~isnan(minDepVinf)
    [iD, jD] = ind2sub(size(vInfDepGrid), idxDep);
    fprintf('Minimum departure v_inf = %.6f km/s\n', minDepVinf);
    fprintf('  Departure: %s\n', datestr(depDates(jD)));
    fprintf('  Arrival  : %s\n', datestr(arrDates(iD)));
    fprintf('  TOF      : %.2f days\n\n', TOFGrid(iD,jD));
else
    iD = [];
    jD = [];
end

if ~isnan(minArrVinf)
    [iA, jA] = ind2sub(size(vInfArrGrid), idxArr);
    fprintf('Minimum arrival v_inf = %.6f km/s\n', minArrVinf);
    fprintf('  Departure: %s\n', datestr(depDates(jA)));
    fprintf('  Arrival  : %s\n', datestr(arrDates(iA)));
    fprintf('  TOF      : %.2f days\n\n', TOFGrid(iA,jA));
else
    iA = [];
    jA = [];
end

% ----------------------------------------------------------
% Best-point structs
% ----------------------------------------------------------
bestC3 = [];
bestDep = [];
bestArr = [];

if ~isempty(iC3)
    bestC3.depDate = depDates(jC3);
    bestC3.arrDate = arrDates(iC3);
    bestC3.label = 'Min C3';
end

if ~isempty(iD)
    bestDep.depDate = depDates(jD);
    bestDep.arrDate = arrDates(iD);
    bestDep.label = 'Min dep v_\infty';
end

if ~isempty(iA)
    bestArr.depDate = depDates(jA);
    bestArr.arrDate = arrDates(iA);
    bestArr.label = 'Min arr v_\infty';
end

% ----------------------------------------------------------
% Heatmap-style porkchops
% ----------------------------------------------------------
tofLevelsHeat = 180:30:360;

figure
astro.plot.porkchopPlot(depDates, arrDates, C3Grid, ...
    'Earth-Mars Porkchop: C3', 'C3 [km^2/s^2]', ...
    TOFGrid, tofLevelsHeat, bestC3);

figure
astro.plot.porkchopPlot(depDates, arrDates, vInfDepGrid, ...
    'Earth-Mars Porkchop: Departure v_{\infty}', ...
    'Departure v_{\infty} [km/s]', ...
    TOFGrid, tofLevelsHeat, bestDep);

figure
astro.plot.porkchopPlot(depDates, arrDates, vInfArrGrid, ...
    'Earth-Mars Porkchop: Arrival v_{\infty}', ...
    'Arrival v_{\infty} [km/s]', ...
    TOFGrid, tofLevelsHeat, bestArr);

% ----------------------------------------------------------
% Classical contour-style porkchops
% Mask high values to avoid ugly contour islands
% ----------------------------------------------------------
C3Plot = C3Grid;
C3Plot(C3Plot > 40) = NaN;

vInfArrPlot = vInfArrGrid;
vInfArrPlot(vInfArrPlot > 6) = NaN;

vInfDepPlot = vInfDepGrid;
vInfDepPlot(vInfDepPlot > 6) = NaN;

c3Levels  = 10:2:40;
tofLevels = 180:30:360;
arrLevels = 2.5:0.2:5.5;
depLevels = 3.0:0.2:6.0;

figure
astro.plot.classicalPorkchopPlot( ...
    depDates, arrDates, C3Plot, c3Levels, ...
    TOFGrid, tofLevels, ...
    'Earth-Mars Classical Porkchop: C3', ...
    'C3 [km^2/s^2]', bestC3);

figure
astro.plot.classicalPorkchopPlot( ...
    depDates, arrDates, vInfDepPlot, depLevels, ...
    TOFGrid, tofLevels, ...
    'Earth-Mars Classical Porkchop: Departure v_{\infty}', ...
    'Departure v_{\infty} [km/s]', bestDep);

figure
astro.plot.classicalPorkchopPlot( ...
    depDates, arrDates, vInfArrPlot, arrLevels, ...
    TOFGrid, tofLevels, ...
    'Earth-Mars Classical Porkchop: Arrival v_{\infty}', ...
    'Arrival v_{\infty} [km/s]', bestArr);