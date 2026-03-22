function sol = solveIzzo(r1, r2, tof, mu, longWay)
%SOLVEIZZO Placeholder interface for Izzo-style Lambert solver.
%
% This function is intentionally structured as the future home of a true
% Izzo Lambert implementation. For now, it uses the already validated
% universal-variable solver internally so that:
%   1) the toolbox API is stabilized,
%   2) examples/tests can be built now,
%   3) the internal algorithm can be upgraded later without changing users'
%      scripts.
%
% INPUTS
%   r1      : initial position vector [km]
%   r2      : final position vector [km]
%   tof     : time of flight [s]
%   mu      : gravitational parameter [km^3/s^2]
%   longWay : logical flag
%             false -> short-way transfer
%             true  -> long-way transfer
%
% OUTPUT
%   sol : struct with fields
%       .v1
%       .v2
%       .converged
%       .iterations
%       .tofError
%       .method
%       .backend
%       .notes

if nargin < 5
    longWay = false;
end

% -------------------------------------------------------------------------
% TEMPORARY BACKEND
% -------------------------------------------------------------------------
% For now, use the validated universal-variable solver internally.
% Later, replace this section with the true Izzo x-lambda formulation.
base = astro.lambert.solveUniversal(r1, r2, tof, mu, longWay);

sol.v1         = base.v1;
sol.v2         = base.v2;
sol.converged  = base.converged;
sol.iterations = base.iterations;
sol.tofError   = base.tofError;

sol.method  = 'Izzo';
sol.backend = 'universal-placeholder';
sol.notes   = ['Temporary Izzo interface using the validated universal ' ...
               'solver backend. Replace internal core with true Izzo ' ...
               'iteration later.'];
end