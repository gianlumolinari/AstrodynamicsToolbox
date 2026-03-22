clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/astroToolbox')
startup

% ==========================================================
% Earth-Mars porkchop plot
% Repo-inspired version:
%   - contour lines only
%   - C3 contours
%   - thick magenta TOF contours
%   - C3 < 31 km^2/s^2 masking
%   - 2005-2007 verification window
% ==========================================================

sun = astro.bodies.getBody('sun');

% ----------------------------------------------------------
% Verification window (same style as reference repo)
% ----------------------------------------------------------
departure_dates = datetime(2005,6,20,0,0,0):days(2):datetime(2005,11,7,0,0,0);
arrival_dates   = datetime(2005,12,1,0,0,0):days(5):datetime(2007,2,24,0,0,0);

numDepartures = numel(departure_dates);
numArrivals   = numel(arrival_dates);

% Set true only if you have Parallel Computing Toolbox
USE_PARFOR = true;

% ----------------------------------------------------------
% Cache Horizons states
% ----------------------------------------------------------
cacheFile = fullfile('data', 'earth_mars_porkchop_cache_2005_2007.mat');
useCache = false;

if exist(cacheFile, 'file')
    S = load(cacheFile);
    if isequal(S.departure_dates, departure_dates) && isequal(S.arrival_dates, arrival_dates)
        earthStates = S.earthStates;
        marsStates  = S.marsStates;
        useCache = true;
        disp('Loaded Horizons states from cache.')
    end
end

if ~useCache
    disp('Fetching Horizons states...')
    
    earthStates = cell(1, numDepartures);
    for i = 1:numDepartures
        epochStr = datestr(departure_dates(i), 'yyyy-mm-dd HH:MM:SS');
        earthStates{i} = astro.ephem.getHorizonsState('399', epochStr, '500@10');
    end
    
    marsStates = cell(1, numArrivals);
    for j = 1:numArrivals
        epochStr = datestr(arrival_dates(j), 'yyyy-mm-dd HH:MM:SS');
        marsStates{j} = astro.ephem.getHorizonsState('499', epochStr, '500@10');
    end
    
    save(cacheFile, 'departure_dates', 'arrival_dates', 'earthStates', 'marsStates');
    disp('Saved Horizons states to cache.')
end

% ----------------------------------------------------------
% Allocate porkchop arrays
% 3rd dimension kept for compatibility with repo structure
% even though we are using one Izzo branch here
% ----------------------------------------------------------
C3  = NaN(numArrivals, numDepartures, 1);
TOF = NaN(numArrivals, numDepartures, 1);

v_inf_1 = NaN(numArrivals, numDepartures, 1);
v_inf_2 = NaN(numArrivals, numDepartures, 1);
vv1     = NaN(numArrivals, numDepartures, 1);
vv2     = NaN(numArrivals, numDepartures, 1);

% ----------------------------------------------------------
% Main Lambert scan
% ----------------------------------------------------------
disp('Running Lambert scan...')

nPairs = numDepartures * numArrivals;

C3_vec      = NaN(nPairs,1);
TOF_vec     = NaN(nPairs,1);
vInf1_vec   = NaN(nPairs,1);
vInf2_vec   = NaN(nPairs,1);
vv1_vec     = NaN(nPairs,1);
vv2_vec     = NaN(nPairs,1);

if USE_PARFOR
    if isempty(gcp('nocreate'))
        parpool;
    end
    
    parfor idx = 1:nPairs
        [j, i] = ind2sub([numArrivals, numDepartures], idx);
        
        depUTC = departure_dates(i);
        arrUTC = arrival_dates(j);
        
        tofSec = seconds(arrUTC - depUTC);
        
        if tofSec <= 0
            continue
        end
        
        D_earth = earthStates{i};
        D_mars  = marsStates{j};
        
        r1 = D_earth.r;
        r2 = D_mars.r;
        
        vE = D_earth.v;
        vM = D_mars.v;
        
        try
            sol = astro.lambert.solveIzzo(r1, r2, tofSec, sun.mu, false);
            if ~sol.converged
                continue
            end
            
            vinf1 = norm(sol.v1 - vE);
            vinf2 = norm(vM - sol.v2);
            
            C3_vec(idx)    = vinf1^2;
            TOF_vec(idx)   = tofSec / 86400;
            vInf1_vec(idx) = vinf1;
            vInf2_vec(idx) = vinf2;
            vv1_vec(idx)   = norm(sol.v1);
            vv2_vec(idx)   = norm(sol.v2);
        catch
        end
    end
else
    for idx = 1:nPairs
        [j, i] = ind2sub([numArrivals, numDepartures], idx);
        
        depUTC = departure_dates(i);
        arrUTC = arrival_dates(j);
        
        tofSec = seconds(arrUTC - depUTC);
        
        if tofSec <= 0
            continue
        end
        
        D_earth = earthStates{i};
        D_mars  = marsStates{j};
        
        r1 = D_earth.r;
        r2 = D_mars.r;
        
        vE = D_earth.v;
        vM = D_mars.v;
        
        try
            sol = astro.lambert.solveIzzo(r1, r2, tofSec, sun.mu, false);
            if ~sol.converged
                continue
            end
            
            vinf1 = norm(sol.v1 - vE);
            vinf2 = norm(vM - sol.v2);
            
            C3_vec(idx)    = vinf1^2;
            TOF_vec(idx)   = tofSec / 86400;
            vInf1_vec(idx) = vinf1;
            vInf2_vec(idx) = vinf2;
            vv1_vec(idx)   = norm(sol.v1);
            vv2_vec(idx)   = norm(sol.v2);
        catch
        end
    end
end

% Reshape into repo-like arrays
C3(:,:,1)      = reshape(C3_vec,    [numArrivals, numDepartures]);
TOF(:,:,1)     = reshape(TOF_vec,   [numArrivals, numDepartures]);
v_inf_1(:,:,1) = reshape(vInf1_vec, [numArrivals, numDepartures]);
v_inf_2(:,:,1) = reshape(vInf2_vec, [numArrivals, numDepartures]);
vv1(:,:,1)     = reshape(vv1_vec,   [numArrivals, numDepartures]);
vv2(:,:,1)     = reshape(vv2_vec,   [numArrivals, numDepartures]);

% ----------------------------------------------------------
% Plot-prep masks
% ----------------------------------------------------------
C3_plot = C3;
C3_plot(C3_plot >= 31) = NaN;

TOF_plot = TOF;
for i = 1:size(C3_plot,1)
    for j = 1:size(C3_plot,2)
        for k = 1:size(C3_plot,3)
            if isnan(C3_plot(i,j,k))
                TOF_plot(i,j,k) = NaN;
            end
        end
    end
end

% ----------------------------------------------------------
% Porkchop plot
% ----------------------------------------------------------
figure('Color','w');
[X, Y] = meshgrid(datenum(departure_dates), datenum(arrival_dates));

levels = 10:2:30;        % C3 [km^2/s^2]
levels_tof = 150:25:600; % TOF [days]

hold on
for k = 1:size(C3,3)
    contour(X, Y, C3_plot(:,:,k), levels, 'LineWidth', 1.2, 'ShowText', 'on');
    contour(X, Y, TOF_plot(:,:,k), levels_tof, 'm-', 'LineWidth', 3, 'ShowText', 'on');
end
hold off

colormap(jet);
cb = colorbar;
ylabel(cb, 'C3 [km^2/s^2]');
clim([10 30]);

ax = gca;
ax.Color = 'w';
ax.XLim = datenum([departure_dates(1), departure_dates(end)]);
ax.YLim = datenum([arrival_dates(1), arrival_dates(end)]);

ax.XTick = linspace(ax.XLim(1), ax.XLim(2), 20);
ax.YTick = linspace(ax.YLim(1), ax.YLim(2), 20);

datetick(ax, 'x', 'yyyy-mm-dd', 'keeplimits', 'keepticks');
datetick(ax, 'y', 'yyyy-mm-dd', 'keeplimits', 'keepticks');

ax.XTickLabelRotation = 45;
ax.FontSize = 11;
ax.FontWeight = 'bold';
ax.XColor = 'k';
ax.YColor = 'k';
ax.LineWidth = 1.2;

set(ax, 'XGrid', 'on', 'YGrid', 'on', ...
    'GridColor', [0.75 0.75 0.75], ...
    'GridAlpha', 0.35);

box on

xlabel('Departure Dates', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k');
ylabel('Arrival Dates', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k');
title('Earth–Mars Pork Chop: C_3 < 30 km^2/s^2', ...
    'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');

% ----------------------------------------------------------
% Print minimum C3 in plotted region
% ----------------------------------------------------------
[minC3, idxMin] = min(C3_plot(:), [], 'omitnan');
if ~isnan(minC3)
    [iMin, jMin, ~] = ind2sub(size(C3_plot), idxMin);
    fprintf('Minimum plotted C3 found = %.6f km^2/s^2\n', minC3);
    fprintf('  Departure: %s\n', datestr(departure_dates(jMin)));
    fprintf('  Arrival  : %s\n', datestr(arrival_dates(iMin)));
    fprintf('  TOF      : %.2f days\n', TOF(iMin,jMin,1));
    fprintf('  Departure v_inf : %.6f km/s\n', v_inf_1(iMin,jMin,1));
    fprintf('  Arrival   v_inf : %.6f km/s\n', v_inf_2(iMin,jMin,1));
end