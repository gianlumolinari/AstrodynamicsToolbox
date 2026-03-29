function out = validateHaloFromLiterature(mu, libr, branch, ref, opts)
%VALIDATEHALOFROMLITERATURE Validate a halo orbit against literature data.
%
% INPUTS
%   mu     : CR3BP mass parameter
%   libr   : 1 or 2
%   branch : 'north' or 'south'
%   ref    : struct with optional fields
%       .Az          target vertical amplitude
%       .Ax          target x-amplitude relative to libration point
%       .period      target period
%       .jacobi      target Jacobi constant
%       .state0      full reference state [6x1]
%       .source      descriptive string
%
%   opts   : optional struct with fields
%       .x0Guess     initial halo x guess
%       .vy0Guess    initial halo vy guess
%       .thalfGuess  initial half-period guess
%       .maxIter     correction iterations
%       .tol         correction tolerance
%       .verbose     true/false
%
% OUTPUT
%   out : struct with fields
%       .ref
%       .corr
%       .orbit
%       .measured
%       .errors
%       .summaryTable

if nargin < 5 || isempty(opts)
    opts = struct();
end
if ~isfield(opts,'maxIter') || isempty(opts.maxIter)
    opts.maxIter = 30;
end
if ~isfield(opts,'tol') || isempty(opts.tol)
    opts.tol = 1e-12;
end
if ~isfield(opts,'verbose') || isempty(opts.verbose)
    opts.verbose = true;
end

branch = lower(string(branch));
if branch ~= "north" && branch ~= "south"
    error('branch must be "north" or "south".');
end

if ~isfield(ref,'source') || isempty(ref.source)
    ref.source = 'literature reference';
end

L = astro.cr3bp.lagrangePoints(mu);
if libr == 1
    xL = L.L1(1);
    x0Default = xL - 1e-2;
elseif libr == 2
    xL = L.L2(1);
    x0Default = xL + 1e-2;
else
    error('Only libr = 1 or 2 supported.');
end

if ~isfield(opts,'x0Guess') || isempty(opts.x0Guess)
    opts.x0Guess = x0Default;
end
if ~isfield(opts,'vy0Guess') || isempty(opts.vy0Guess)
    opts.vy0Guess = 0.1;
end
if ~isfield(opts,'thalfGuess') || isempty(opts.thalfGuess)
    opts.thalfGuess = 1.5;
end

% -------------------------------------------------------------------------
% Build corrected orbit
% -------------------------------------------------------------------------
if isfield(ref,'state0') && ~isempty(ref.state0)
    stateRef = ref.state0(:);

    corr = astro.cr3bp.differentialCorrectionHalo( ...
        stateRef(1), stateRef(3), stateRef(5), 0.5*ref.period, ...
        mu, opts.maxIter, opts.tol);

else
    if ~isfield(ref,'Az') || isempty(ref.Az)
        error('For halo literature validation, provide either ref.state0 or ref.Az.');
    end

    if branch == "north"
        z0Guess = +abs(ref.Az);
    else
        z0Guess = -abs(ref.Az);
    end

    corr = astro.cr3bp.differentialCorrectionHalo( ...
        opts.x0Guess, z0Guess, opts.vy0Guess, opts.thalfGuess, ...
        mu, opts.maxIter, opts.tol);
end

if ~corr.converged
    error('Halo correction did not converge.');
end

orbit = astro.periodic.buildOrbitStruct(corr.state0, corr.period, mu, ...
    struct('family','halo', ...
           'libr', libr, ...
           'system', 'CR3BP', ...
           'source', ref.source, ...
           'dimension', '3D'), ...
    struct('verbose', opts.verbose));

x = orbit.x(:,1);
z = orbit.x(:,3);

measured = struct();
measured.x0 = corr.state0(1);
measured.z0 = corr.state0(3);
measured.vy0 = corr.state0(5);
measured.period = corr.period;
measured.jacobi = astro.cr3bp.jacobiConstant(corr.state0.', mu);
measured.Ax = max(abs(x - xL));
measured.Az = max(abs(z));
measured.closureError = orbit.closureError;

errors = struct();
errors.Ax = localDiff(ref, measured, 'Ax');
errors.Az = localDiff(ref, measured, 'Az');
errors.period = localDiff(ref, measured, 'period');
errors.jacobi = localDiff(ref, measured, 'jacobi');

out = struct();
out.ref = ref;
out.corr = corr;
out.orbit = orbit;
out.measured = measured;
out.errors = errors;

names = {};
targets = [];
values = [];
absErrs = [];

[names, targets, values, absErrs] = localAppend(names, targets, values, absErrs, ref, measured, errors, 'Ax');
[names, targets, values, absErrs] = localAppend(names, targets, values, absErrs, ref, measured, errors, 'Az');
[names, targets, values, absErrs] = localAppend(names, targets, values, absErrs, ref, measured, errors, 'period');
[names, targets, values, absErrs] = localAppend(names, targets, values, absErrs, ref, measured, errors, 'jacobi');

if isempty(names)
    out.summaryTable = table();
else
    out.summaryTable = table(string(names(:)), targets(:), values(:), absErrs(:), ...
        'VariableNames', {'Quantity','Target','Measured','AbsError'});
end

if opts.verbose
    fprintf('\nHalo literature validation\n');
    fprintf('--------------------------\n');
    fprintf('Source           : %s\n', ref.source);
    fprintf('Libration point  : L%d\n', libr);
    fprintf('Branch           : %s\n', branch);
    fprintf('xL               : %.16f\n', xL);
    fprintf('x0 corrected     : %.16f\n', measured.x0);
    fprintf('z0 corrected     : %.16f\n', measured.z0);
    fprintf('vy0 corrected    : %.16f\n', measured.vy0);
    fprintf('Period measured  : %.16f\n', measured.period);
    fprintf('Jacobi measured  : %.16f\n', measured.jacobi);
    fprintf('Ax measured      : %.16e\n', measured.Ax);
    fprintf('Az measured      : %.16e\n', measured.Az);
    fprintf('Closure error    : %.3e\n', measured.closureError);

    if ~isempty(out.summaryTable)
        disp(out.summaryTable)
    end
end
end

function err = localDiff(ref, measured, fieldName)
if isfield(ref, fieldName) && ~isempty(ref.(fieldName))
    err = abs(measured.(fieldName) - ref.(fieldName));
else
    err = NaN;
end
end

function [names, targets, values, absErrs] = localAppend(names, targets, values, absErrs, ref, measured, errors, fieldName)
if isfield(ref, fieldName) && ~isempty(ref.(fieldName))
    names{end+1} = fieldName; %#ok<AGROW>
    targets(end+1) = ref.(fieldName); %#ok<AGROW>
    values(end+1) = measured.(fieldName); %#ok<AGROW>
    absErrs(end+1) = errors.(fieldName); %#ok<AGROW>
end
end