function out = generateHaloByAmplitude(mu, libr, AzTarget, branch, opts)
%GENERATEHALOBYAMPLITUDE Generate a halo orbit from a target vertical amplitude.
%
% INPUTS
%   mu       : CR3BP mass parameter
%   libr     : 1 or 2
%   AzTarget : target vertical amplitude seed, used as z0 initial design parameter
%   branch   : 'north' or 'south'
%   opts     : optional struct with fields
%       .x0Guess        : optional initial guess for x0
%       .vy0Guess       : optional initial guess for vy0
%       .thalfGuess     : optional initial guess for half period
%       .maxIter        : correction iterations, default 30
%       .tol            : correction tolerance, default 1e-12
%       .buildOrbit     : true/false, default true
%       .verbose        : true/false, default true
%
% OUTPUT
%   out : struct with fields
%       .libr
%       .branch
%       .xL
%       .AzTarget
%       .z0Guess
%       .x0Guess
%       .vy0Guess
%       .thalfGuess
%       .corr
%       .state0
%       .period
%       .orbit
%       .AzMeasured
%       .AxMeasured

if nargin < 5 || isempty(opts)
    opts = struct();
end

branch = lower(string(branch));
if branch ~= "north" && branch ~= "south"
    error('branch must be "north" or "south".');
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

if ~isfield(opts, 'x0Guess') || isempty(opts.x0Guess)
    if libr == 1
        opts.x0Guess = xL - 1e-2;
    else
        opts.x0Guess = xL + 1e-2;
    end
end
if ~isfield(opts, 'vy0Guess') || isempty(opts.vy0Guess)
    opts.vy0Guess = 0.1;
end
if ~isfield(opts, 'thalfGuess') || isempty(opts.thalfGuess)
    opts.thalfGuess = 1.5;
end

if branch == "north"
    z0Guess = +abs(AzTarget);
else
    z0Guess = -abs(AzTarget);
end

if opts.verbose
    fprintf('\nGenerating halo orbit by amplitude\n');
    fprintf('---------------------------------\n');
    fprintf('Libration point : L%d\n', libr);
    fprintf('Branch          : %s\n', branch);
    fprintf('xL              : %.16f\n', xL);
    fprintf('AzTarget        : %.16e\n', AzTarget);
    fprintf('x0Guess         : %.16f\n', opts.x0Guess);
    fprintf('z0Guess         : %.16f\n', z0Guess);
    fprintf('vy0Guess        : %.16f\n', opts.vy0Guess);
    fprintf('thalfGuess      : %.16f\n', opts.thalfGuess);
end

corr = astro.cr3bp.differentialCorrectionHalo( ...
    opts.x0Guess, z0Guess, opts.vy0Guess, opts.thalfGuess, ...
    mu, opts.maxIter, opts.tol);

if ~corr.converged
    error('Halo correction did not converge.');
end

out = struct();
out.libr = libr;
out.branch = char(branch);
out.xL = xL;
out.AzTarget = AzTarget;
out.z0Guess = z0Guess;
out.x0Guess = opts.x0Guess;
out.vy0Guess = opts.vy0Guess;
out.thalfGuess = opts.thalfGuess;
out.corr = corr;
out.state0 = corr.state0;
out.period = corr.period;
out.orbit = [];
out.AzMeasured = max(abs(corr.state0(3)));
out.AxMeasured = abs(corr.state0(1) - xL);

if opts.buildOrbit
    orbit = astro.periodic.buildOrbitStruct(corr.state0, corr.period, mu, ...
        struct('family', 'halo', ...
               'libr', libr, ...
               'system', 'CR3BP', ...
               'source', 'generateHaloByAmplitude', ...
               'dimension', '3D'), ...
        struct('verbose', opts.verbose));

    x = orbit.x(:,1);
    z = orbit.x(:,3);

    out.orbit = orbit;
    out.AzMeasured = max(abs(z));
    out.AxMeasured = max(abs(x - xL));

    if opts.verbose
        fprintf('\nCorrected halo orbit\n');
        fprintf('--------------------\n');
        fprintf('x0              : %.16f\n', corr.state0(1));
        fprintf('z0              : %.16f\n', corr.state0(3));
        fprintf('vy0             : %.16f\n', corr.state0(5));
        fprintf('T               : %.16f\n', corr.period);
        fprintf('closure error    : %.3e\n', orbit.closureError);
        fprintf('Az measured      : %.16e\n', out.AzMeasured);
        fprintf('Ax measured      : %.16e\n', out.AxMeasured);
    end
end
end