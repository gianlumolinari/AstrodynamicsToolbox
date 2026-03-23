function result = propagateMGA(sequence, dates, muSun, backend)
%PROPAGATEMGA Simple multi-gravity-assist leg evaluation.
%
% INPUTS
%   sequence : cell array of body names, e.g.
%              {'earth','venus','earth','jupiter'}
%   dates    : cell array of UTC date strings, one per body encounter
%   muSun    : Sun gravitational parameter [km^3/s^2]
%   backend  : 'spice' or 'horizons' (default: 'spice')
%
% OUTPUT
%   result : struct containing leg-by-leg Lambert information
%
% NOTES
%   This first version:
%   - solves each heliocentric leg independently
%   - computes departure/arrival v_infinity for each encounter
%   - does not yet enforce flyby continuity constraints automatically

if nargin < 4 || isempty(backend)
    backend = 'spice';
end

if numel(sequence) ~= numel(dates)
    error('sequence and dates must have the same length.');
end

nBodies = numel(sequence);
nLegs = nBodies - 1;

result.sequence = sequence;
result.dates = dates;
result.backend = backend;
result.legs = struct([]);

for k = 1:nLegs
    body1 = lower(string(sequence{k}));
    body2 = lower(string(sequence{k+1}));

    epoch1 = char(string(dates{k}));
    epoch2 = char(string(dates{k+1}));

    t1 = datetime(epoch1, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
    t2 = datetime(epoch2, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
    tofSec = seconds(t2 - t1);

    if tofSec <= 0
        error('Each leg must have strictly positive time of flight.');
    end

    state1 = astro.ephem.getState(backend, body1, epoch1);
    state2 = astro.ephem.getState(backend, body2, epoch2);

    sol = astro.lambert.solveIzzo(state1.r, state2.r, tofSec, muSun, false);

    if ~sol.converged
        error('Lambert solver failed on leg %d: %s -> %s', k, body1, body2);
    end

    vInfDep = sol.v1 - state1.v;
    vInfArr = state2.v - sol.v2;

    result.legs(k).bodyDepart = char(body1);
    result.legs(k).bodyArrive = char(body2);
    result.legs(k).epochDepart = epoch1;
    result.legs(k).epochArrive = epoch2;
    result.legs(k).tofSec = tofSec;
    result.legs(k).tofDays = tofSec / 86400;

    result.legs(k).r1 = state1.r;
    result.legs(k).r2 = state2.r;
    result.legs(k).vPlanet1 = state1.v;
    result.legs(k).vPlanet2 = state2.v;

    result.legs(k).v1 = sol.v1;
    result.legs(k).v2 = sol.v2;

    result.legs(k).vInfDep = vInfDep;
    result.legs(k).vInfArr = vInfArr;

    result.legs(k).vInfDepMag = norm(vInfDep);
    result.legs(k).vInfArrMag = norm(vInfArr);

    result.legs(k).lambert = sol;
end
end