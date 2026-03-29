function out = generatePlanarLyapunovByAmplitude(mu, libr, AxTarget, opts)
%GENERATEPLANARLYAPUNOVBYAMPLITUDE Generate a planar Lyapunov orbit from a target in-plane amplitude.
%
% INPUTS
%   mu       : CR3BP mass parameter
%   libr     : 1 or 2
%   AxTarget : target in-plane amplitude relative to the chosen libration point
%              AxTarget = |x0 - xL|
%   opts     : optional struct with fields
%       .xSign          : sign convention for x0 offset; default chosen automatically
%       .vy0Guess       : optional initial guess for vy0
%       .maxIter        : correction iterations, default 30
%       .tol            : correction tolerance, default 1e-12
%       .buildOrbit     : true/false, default true
%       .verbose        : true/false, default true
%
% OUTPUT
%   out : struct with fields
%       .libr
%       .xL
%       .AxTarget
%       .x0Guess
%       .vy0Guess
%       .corr
%       .state0
%       .period
%       .orbit
%       .AxMeasured
%       .AyMeasured

if nargin < 4 || isempty(opts)
    opts = struct();
end

if ~isfield(opts, 'xSign') || isempty(opts.xSign)
    if libr == 1
        opts.xSign = -1;
    elseif libr == 2
        opts.xSign = +1;
    else
        error('Only libr = 1 or 2 supported.');
    end
end

if ~isfield(opts, 'vy0Guess') || isempty(opts.vy0Guess)
    opts.vy0Guess = 0.1;
end
if ~isfield(opts, 'maxIter') || isempty(opts.maxIter)
    opts.maxIter = 30;
end
if ~isfield(opts, 'tol') || isempty(opts.tol)
    opts.tol = 1e-12;
end
if ~isfield(opts, 'buildOrbit') || isempty(opts.buildOrbit)
    opts.buildOrbit = true;
end
if ~isfield(opts, 'verbose') || isempty(opts.verbose)
    opts.verbose = true;
end

L = astro.cr3bp.lagrangePoints(mu);
if libr == 1
    xL = L.L1(1);
elseif libr == 2
    xL = L.L2(1);
else
    error('Only libr = 1 or 2 supported.');
end

x0Guess = xL + opts.xSign * AxTarget;
vy0Guess = opts.vy0Guess;

if opts.verbose
    fprintf('\nGenerating planar Lyapunov orbit by amplitude\n');
    fprintf('--------------------------------------------\n');
    fprintf('Libration point : L%d\n', libr);
    fprintf('xL              : %.16f\n', xL);
    fprintf('AxTarget        : %.16e\n', AxTarget);
    fprintf('x0Guess         : %.16f\n', x0Guess);
    fprintf('vy0Guess        : %.16f\n', vy0Guess);
end

corr = astro.cr3bp.differentialCorrectionPlanarLyapunov( ...
    x0Guess, vy0Guess, mu, opts.maxIter, opts.tol);

if ~corr.converged
    error('Planar Lyapunov correction did not converge.');
end

out = struct();
out.libr = libr;
out.xL = xL;
out.AxTarget = AxTarget;
out.x0Guess = x0Guess;
out.vy0Guess = vy0Guess;
out.corr = corr;
out.state0 = corr.state0;
out.period = corr.period;
out.orbit = [];
out.AxMeasured = abs(corr.state0(1) - xL);
out.AyMeasured = NaN;

if opts.buildOrbit
    orbit = astro.periodic.buildOrbitStruct(corr.state0, corr.period, mu, ...
        struct('family', 'planar_lyapunov', ...
               'libr', libr, ...
               'system', 'CR3BP', ...
               'source', 'generatePlanarLyapunovByAmplitude', ...
               'dimension', '2D'), ...
        struct('verbose', opts.verbose));

    x = orbit.x(:,1);
    y = orbit.x(:,2);

    out.orbit = orbit;
    out.AxMeasured = max(abs(x - xL));
    out.AyMeasured = max(abs(y));

    if opts.verbose
        fprintf('\nCorrected planar Lyapunov orbit\n');
        fprintf('-------------------------------\n');
        fprintf('x0              : %.16f\n', corr.state0(1));
        fprintf('vy0             : %.16f\n', corr.state0(5));
        fprintf('T               : %.16f\n', corr.period);
        fprintf('closure error    : %.3e\n', orbit.closureError);
        fprintf('Ax measured      : %.16e\n', out.AxMeasured);
        fprintf('Ay measured      : %.16e\n', out.AyMeasured);
    end
end
end