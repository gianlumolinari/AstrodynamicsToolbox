function results = scanMGAWindow(sequence, depDatesUTC, tofGridDays, muSun, bodyMap, minAltMap, backend)
%SCANMGAWINDOW Coarse MGA window scan over departure dates and leg TOFs.
%
% INPUTS
%   sequence    : cell array of bodies
%   depDatesUTC : cell array of departure UTC strings
%   tofGridDays : cell array, one vector per leg
%   muSun       : Sun gravitational parameter
%   bodyMap     : struct of flyby body structs
%   minAltMap   : struct of minimum flyby altitudes [km]
%   backend     : 'spice' or 'horizons'
%
% OUTPUT
%   results : struct array with fields
%       .depUTC
%       .tofDays
%       .dates
%       .objective
%       .maxMagMismatch
%       .maxTurnViolationDeg
%       .allFlybysFeasible
%
% NOTES
%   This is a brute-force coarse scanner meant to identify promising
%   windows before local optimization.

if nargin < 7 || isempty(backend)
    backend = 'spice';
end

nLegs = numel(sequence) - 1;
if numel(tofGridDays) ~= nLegs
    error('tofGridDays must contain one grid per leg.');
end

% Build all TOF combinations
gridCells = tofGridDays(:).';
[gridMesh{1:nLegs}] = ndgrid(gridCells{:});
nComb = numel(gridMesh{1});

tofComb = zeros(nComb, nLegs);
for k = 1:nLegs
    tofComb(:,k) = gridMesh{k}(:);
end

results = struct([]);
idxRes = 0;

for iDep = 1:numel(depDatesUTC)
    depUTC = depDatesUTC{iDep};

    for iComb = 1:nComb
        tofDays = tofComb(iComb,:).';

        try
            dates = astro.mga.unpackDatesFromTOF(depUTC, tofDays);
            traj = astro.mga.propagateMGA(sequence, dates, muSun, backend);
            report = astro.mga.evaluateMGA(traj, bodyMap, minAltMap, 1e-3);
            J = astro.mga.mgaObjective(tofDays, sequence, depUTC, muSun, bodyMap, minAltMap, backend);

            maxMagMismatch = 0;
            maxTurnViolationDeg = 0;
            allFlybysFeasible = true;

            for k = 1:report.nFlybys
                chk = report.flybys(k).check;
                maxMagMismatch = max(maxMagMismatch, chk.magMismatch);
                maxTurnViolationDeg = max(maxTurnViolationDeg, chk.deltaReqDeg - chk.deltaMaxDeg);
                allFlybysFeasible = allFlybysFeasible && chk.feasible;
            end

            idxRes = idxRes + 1;
            results(idxRes).depUTC = depUTC;
            results(idxRes).tofDays = tofDays;
            results(idxRes).dates = dates;
            results(idxRes).objective = J;
            results(idxRes).maxMagMismatch = maxMagMismatch;
            results(idxRes).maxTurnViolationDeg = maxTurnViolationDeg;
            results(idxRes).allFlybysFeasible = allFlybysFeasible;

        catch
        end
    end
end

if ~isempty(results)
    [~, I] = sort([results.objective]);
    results = results(I);
end
end