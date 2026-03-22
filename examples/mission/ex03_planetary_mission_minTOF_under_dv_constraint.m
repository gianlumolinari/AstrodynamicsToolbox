clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/astroToolbox')
startup

% ==========================================================
% Generic planetary mission search
% Objective:
%   minimize time of flight subject to total delta-v limit
%
% Primary body sequence:
%   Earth departure -> target planet arrival
%
% Solver:
%   Izzo Lambert
% ==========================================================

sun   = astro.bodies.getBody('sun');
earth = astro.bodies.getBody('earth');

% ----------------------------------------------------------
% USER SETTINGS
% ----------------------------------------------------------
targetName = 'venus';
dvMax = 10.0;

% Parking orbit assumptions
hDep = 300;              % km, Earth parking orbit altitude
hArr = 300;              % km, target capture orbit altitude

rDep = earth.radius + hDep;

% Date windows
departure_dates = datetime(2026,9,1,0,0,0):days(2):datetime(2026,12,1,0,0,0);
arrival_dates   = datetime(2027,1,1,0,0,0):days(3):datetime(2028,6,1,0,0,0);

% Minimum time of flight to avoid nonsense cases
minTOFdays = 50;

% Optional parallel scan
USE_PARFOR = true;

% ----------------------------------------------------------
% TARGET BODY SETUP
% ----------------------------------------------------------
target = astro.bodies.getBody(targetName);
rArr = target.radius + hArr;

switch lower(targetName)
    case 'mercury'
        targetID = '199';
    case 'venus'
        targetID = '299';
    case 'earth'
        targetID = '399';
    case 'mars'
        targetID = '499';
    case 'jupiter'
        targetID = '599';
    case 'saturn'
        targetID = '699';
    case 'uranus'
        targetID = '799';
    case 'neptune'
        targetID = '899';
    otherwise
        error('Unsupported targetName: %s', targetName);
end

numDepartures = numel(departure_dates);
numArrivals   = numel(arrival_dates);

% ----------------------------------------------------------
% Cache Horizons states
% ----------------------------------------------------------
cacheFile = fullfile('data', ['mission_cache_earth_to_' lower(targetName) '.mat']);
useCache = false;

if exist(cacheFile, 'file')
    S = load(cacheFile);
    if isequal(S.departure_dates, departure_dates) && isequal(S.arrival_dates, arrival_dates)
        earthStates  = S.earthStates;
        targetStates = S.targetStates;
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
    
    targetStates = cell(1, numArrivals);
    for j = 1:numArrivals
        epochStr = datestr(arrival_dates(j), 'yyyy-mm-dd HH:MM:SS');
        targetStates{j} = astro.ephem.getHorizonsState(targetID, epochStr, '500@10');
    end
    
    save(cacheFile, 'departure_dates', 'arrival_dates', 'earthStates', 'targetStates');
    disp('Saved Horizons states to cache.')
end

% ----------------------------------------------------------
% Allocate storage
% ----------------------------------------------------------
nPairs = numDepartures * numArrivals;

TOF_vec    = NaN(nPairs,1);
C3_vec     = NaN(nPairs,1);
vInfDep_vec = NaN(nPairs,1);
vInfArr_vec = NaN(nPairs,1);
depDV_vec  = NaN(nPairs,1);
arrDV_vec  = NaN(nPairs,1);
totDV_vec  = NaN(nPairs,1);

disp('Running constrained mission scan...')

% ----------------------------------------------------------
% Main scan
% ----------------------------------------------------------
if USE_PARFOR
    if isempty(gcp('nocreate'))
        parpool;
    end
    
    parfor idx = 1:nPairs
        [j, i] = ind2sub([numArrivals, numDepartures], idx);
        
        depUTC = departure_dates(i);
        arrUTC = arrival_dates(j);
        
        tofSec = seconds(arrUTC - depUTC);
        tofDays = tofSec / 86400;
        
        if tofDays <= minTOFdays
            continue
        end
        
        D_earth = earthStates{i};
        D_target = targetStates{j};
        
        r1 = D_earth.r;
        r2 = D_target.r;
        vE = D_earth.v;
        vT = D_target.v;
        
        try
            sol = astro.lambert.solveIzzo(r1, r2, tofSec, sun.mu, false);
            if ~sol.converged
                continue
            end
            
            vInfDep = norm(sol.v1 - vE);
            vInfArr = norm(vT - sol.v2);
            
            dep = astro.maneuvers.departureDeltaV(earth.mu, rDep, vInfDep);
            arr = astro.maneuvers.arrivalDeltaV(target.mu, rArr, vInfArr);
            
            depDV = dep.deltaV;
            arrDV = arr.deltaV;
            totDV = depDV + arrDV;
            
            TOF_vec(idx)     = tofDays;
            C3_vec(idx)      = vInfDep^2;
            vInfDep_vec(idx) = vInfDep;
            vInfArr_vec(idx) = vInfArr;
            depDV_vec(idx)   = depDV;
            arrDV_vec(idx)   = arrDV;
            totDV_vec(idx)   = totDV;
        catch
        end
    end
else
    for idx = 1:nPairs
        [j, i] = ind2sub([numArrivals, numDepartures], idx);
        
        depUTC = departure_dates(i);
        arrUTC = arrival_dates(j);
        
        tofSec = seconds(arrUTC - depUTC);
        tofDays = tofSec / 86400;
        
        if tofDays <= minTOFdays
            continue
        end
        
        D_earth = earthStates{i};
        D_target = targetStates{j};
        
        r1 = D_earth.r;
        r2 = D_target.r;
        vE = D_earth.v;
        vT = D_target.v;
        
        try
            sol = astro.lambert.solveIzzo(r1, r2, tofSec, sun.mu, false);
            if ~sol.converged
                continue
            end
            
            vInfDep = norm(sol.v1 - vE);
            vInfArr = norm(vT - sol.v2);
            
            dep = astro.maneuvers.departureDeltaV(earth.mu, rDep, vInfDep);
            arr = astro.maneuvers.arrivalDeltaV(target.mu, rArr, vInfArr);
            
            depDV = dep.deltaV;
            arrDV = arr.deltaV;
            totDV = depDV + arrDV;
            
            TOF_vec(idx)     = tofDays;
            C3_vec(idx)      = vInfDep^2;
            vInfDep_vec(idx) = vInfDep;
            vInfArr_vec(idx) = vInfArr;
            depDV_vec(idx)   = depDV;
            arrDV_vec(idx)   = arrDV;
            totDV_vec(idx)   = totDV;
        catch
        end
    end
end

% ----------------------------------------------------------
% Reshape into grids
% ----------------------------------------------------------
TOF     = reshape(TOF_vec,     [numArrivals, numDepartures]);
C3      = reshape(C3_vec,      [numArrivals, numDepartures]);
vInfDep = reshape(vInfDep_vec, [numArrivals, numDepartures]);
vInfArr = reshape(vInfArr_vec, [numArrivals, numDepartures]);
depDV   = reshape(depDV_vec,   [numArrivals, numDepartures]);
arrDV   = reshape(arrDV_vec,   [numArrivals, numDepartures]);
totDV   = reshape(totDV_vec,   [numArrivals, numDepartures]);

% ----------------------------------------------------------
% Feasible region
% ----------------------------------------------------------
feasibleMask = totDV <= dvMax;

TOF_feasible = TOF;
TOF_feasible(~feasibleMask) = NaN;

% ----------------------------------------------------------
% Find optimal feasible mission
% ----------------------------------------------------------
[minTOF, idxOpt] = min(TOF_feasible(:), [], 'omitnan');

if isnan(minTOF)
    fprintf('\nNo feasible Earth-to-%s solution found for dvMax = %.3f km/s\n', ...
        targetName, dvMax);
else
    [jOpt, iOpt] = ind2sub([numArrivals, numDepartures], idxOpt);
    
    fprintf('\nOptimal Earth-to-%s mission under delta-v constraint\n', upper(targetName));
    fprintf('------------------------------------------------------\n');
    fprintf('Delta-v limit              : %.6f km/s\n', dvMax);
    fprintf('Departure date             : %s\n', datestr(departure_dates(iOpt)));
    fprintf('Arrival date               : %s\n', datestr(arrival_dates(jOpt)));
    fprintf('Minimum feasible TOF       : %.2f days\n', TOF(jOpt, iOpt));
    fprintf('Departure C3               : %.6f km^2/s^2\n', C3(jOpt, iOpt));
    fprintf('Departure v_inf            : %.6f km/s\n', vInfDep(jOpt, iOpt));
    fprintf('Arrival v_inf              : %.6f km/s\n', vInfArr(jOpt, iOpt));
    fprintf('Departure delta-v          : %.6f km/s\n', depDV(jOpt, iOpt));
    fprintf('Arrival delta-v            : %.6f km/s\n', arrDV(jOpt, iOpt));
    fprintf('Total delta-v              : %.6f km/s\n', totDV(jOpt, iOpt));
end