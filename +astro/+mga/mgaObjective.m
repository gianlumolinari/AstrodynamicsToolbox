function J = mgaObjective(tofDays, sequence, t0UTC, muSun, bodyMap, minAltMap, backend)
%MGAOBJECTIVE Objective function for simple MGA date optimization.
%
% INPUTS
%   tofDays   : leg times of flight [days]
%   sequence  : cell array of body names
%   t0UTC     : departure UTC string
%   muSun     : Sun gravitational parameter [km^3/s^2]
%   bodyMap   : struct of flyby body structs
%   minAltMap : struct of minimum flyby altitudes [km]
%   backend   : 'spice' or 'horizons'
%
% OUTPUT
%   J         : scalar penalty to minimize
%
% OBJECTIVE STRUCTURE
%   J = sum over flybys of:
%       magnitude mismatch penalty
%     + turning infeasibility penalty
%     + optional launch energy penalty
%
% NOTES
%   This first version is intentionally simple and robust.

if nargin < 7 || isempty(backend)
    backend = 'spice';
end

% Basic guards
if any(tofDays <= 0)
    J = 1e12;
    return
end

try
    dates = astro.mga.unpackDatesFromTOF(t0UTC, tofDays);
    traj = astro.mga.propagateMGA(sequence, dates, muSun, backend);
    report = astro.mga.evaluateMGA(traj, bodyMap, minAltMap, 1e-3);

    J = 0;

    % Launch penalty: discourage absurd launch energies
    launchPenalty = traj.legs(1).vInfDepMag^2;
    J = J + 0.02 * launchPenalty;

    % Flyby continuity / feasibility penalty
    for k = 1:report.nFlybys
        chk = report.flybys(k).check;

        magPenalty = chk.magMismatch^2;

        turnExcessDeg = max(0, chk.deltaReqDeg - chk.deltaMaxDeg);
        turnPenalty = 10 * turnExcessDeg^2;

        J = J + magPenalty + turnPenalty;
    end

    % Mild TOF regularization
    J = J + 1e-4 * sum(tofDays.^2);

catch
    J = 1e12;
end
end