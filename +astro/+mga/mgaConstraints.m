function [c, ceq] = mgaConstraints(tofDays, sequence, t0UTC, muSun, bodyMap, minAltMap, backend)
%MGACONSTRAINTS Nonlinear constraints for ballistic MGA optimization.
%
% Inequality constraints c(x) <= 0:
%   delta_req - delta_max <= 0
%
% Equality constraints ceq(x) = 0:
%   |v_inf^-| - |v_inf^+| = 0
%
% This first version returns one equality and one inequality per flyby.

if nargin < 7 || isempty(backend)
    backend = 'spice';
end

% Default: infeasible
c = 1e6;
ceq = 1e6;

if any(tofDays <= 0)
    return
end

try
    dates = astro.mga.unpackDatesFromTOF(t0UTC, tofDays);
    traj = astro.mga.propagateMGA(sequence, dates, muSun, backend);
    report = astro.mga.evaluateMGA(traj, bodyMap, minAltMap, 1e-3);

    nFlybys = report.nFlybys;
    c = zeros(nFlybys, 1);
    ceq = zeros(nFlybys, 1);

    for k = 1:nFlybys
        chk = report.flybys(k).check;

        % hard equality on v_inf magnitude continuity
        ceq(k) = chk.magIn - chk.magOutReq;

        % hard inequality on turning capability
        c(k) = chk.deltaReqDeg - chk.deltaMaxDeg;
    end

catch
    c = 1e6;
    ceq = 1e6;
end
end