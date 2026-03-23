function report = evaluateMGA(traj, bodyMap, minAltMap, magTol)
%EVALUATEMGA Evaluates flyby continuity and feasibility for an MGA sequence.
%
% INPUTS
%   traj      : output from astro.mga.propagateMGA
%   bodyMap   : struct mapping body names to body structs
%               e.g. bodyMap.venus = astro.bodies.getBody('venus');
%   minAltMap : struct mapping flyby bodies to minimum altitude [km]
%               e.g. minAltMap.venus = 300;
%   magTol    : allowable |v_inf| mismatch [km/s] (default: 1e-3)
%
% OUTPUT
%   report : struct with fields
%       .nFlybys
%       .flybys(k)
%
% Each flyby entry contains:
%   .body
%   .legIn
%   .legOut
%   .vInfIn
%   .vInfOutReq
%   .rpMin
%   .minAltitude
%   .check    (output of flybyFeasibility)

if nargin < 4 || isempty(magTol)
    magTol = 1e-3;
end

nLegs = numel(traj.legs);
nFlybys = max(nLegs - 1, 0);

report.nFlybys = nFlybys;
report.flybys = struct([]);

for k = 1:nFlybys
    bodyName = lower(string(traj.sequence{k+1}));
    bodyField = matlab.lang.makeValidName(char(bodyName));

    if ~isfield(bodyMap, bodyField)
        error('Missing body struct for flyby body: %s', bodyName);
    end
    if ~isfield(minAltMap, bodyField)
        error('Missing minimum altitude for flyby body: %s', bodyName);
    end

    body = bodyMap.(bodyField);
    minAlt = minAltMap.(bodyField);
    rpMin = body.radius + minAlt;

    vInfIn = traj.legs(k).vInfArr;
    vInfOutReq = traj.legs(k+1).vInfDep;

    check = astro.mga.flybyFeasibility(vInfIn, vInfOutReq, body.mu, rpMin, magTol);

    report.flybys(k).body = char(bodyName);
    report.flybys(k).legIn = k;
    report.flybys(k).legOut = k + 1;
    report.flybys(k).vInfIn = vInfIn;
    report.flybys(k).vInfOutReq = vInfOutReq;
    report.flybys(k).rpMin = rpMin;
    report.flybys(k).minAltitude = minAlt;
    report.flybys(k).check = check;
end
end