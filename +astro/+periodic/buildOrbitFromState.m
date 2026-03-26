function orbit = buildOrbitFromState(state0, T, mu, opts)
%BUILDORBITFROMSTATE Propagate one full periodic orbit and package output.
%
% INPUTS
%   state0 : 6x1 initial state
%   T      : period
%   mu     : CR3BP mass parameter
%   opts   : optional struct
%
% OUTPUT
%   orbit : struct with fields
%       .x
%       .t
%       .T
%       .state0
%       .monodromy
%       .closureError

    if nargin < 4 || isempty(opts)
        opts = struct();
    end

    if ~isfield(opts,'RelTol'), opts.RelTol = 1e-12; end
    if ~isfield(opts,'AbsTol'), opts.AbsTol = 1e-12; end
    if ~isfield(opts,'verbose'), opts.verbose = true; end

    prop = astro.cr3bp.propagateWithSTM( ...
        state0, [0 T], mu, odeset('RelTol',opts.RelTol,'AbsTol',opts.AbsTol));

    orbit = struct();
    orbit.x = prop.x(:,1:6);
    orbit.t = prop.t(:);
    orbit.T = T;
    orbit.state0 = state0;
    orbit.monodromy = prop.PhiFinal;

    closureErr = norm(orbit.x(end,:).' - state0);
    orbit.closureError = closureErr;

    if opts.verbose
        fprintf('Orbit closure error after one period: %.3e\n', closureErr);
    end
end