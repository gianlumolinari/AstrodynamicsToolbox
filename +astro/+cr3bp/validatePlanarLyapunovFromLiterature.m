function out = validatePlanarLyapunovFromLiterature(mu, libr, ref, opts)
%VALIDATEPLANARLYAPUNOVFROMLITERATURE Validate a planar Lyapunov orbit against literature data.
%
% INPUTS
%   mu   : CR3BP mass parameter
%   libr : 1 or 2
%   ref  : struct with optional fields
%       .Ax          target x-amplitude relative to libration point
%       .Ay          target max |y|
%       .period      target period
%       .jacobi      target Jacobi constant
%       .state0      full reference state [6x1]
%       .source      descriptive string
%
%   opts : optional struct with fields
%       .vy0Guess    initial guess for planar corrector
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

if nargin < 4 || isempty(opts)
    opts = struct();
end
if ~isfield(opts,'vy0Guess') || isempty(opts.vy0Guess)
    opts.vy0Guess = 0.1;
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

if ~isfield(ref,'source') || isempty(ref.source)
    ref.source = 'literature reference';
end

L = astro.cr3bp.lagrangePoints(mu);
if libr == 1
    xL = L.L1(1);
    xSignDefault = -1;
elseif libr == 2
    xL = L.L2(1);
    xSignDefault = +1;
else
    error('Only libr = 1 or 2 supported.');
end

% -------------------------------------------------------------------------
% Build corrected orbit
% -------------------------------------------------------------------------
if isfield(ref,'state0') && ~isempty(ref.state0)
    stateRef = ref.state0(:);

    corr = astro.cr3bp.differentialCorrectionPlanarLyapunov( ...
        stateRef(1), stateRef(5), mu, opts.maxIter, opts.tol);

else
    if ~isfield(ref,'Ax') || isempty(ref.Ax)
        error('For planar literature validation, provide either ref.state0 or ref.Ax.');
    end

    x0Guess = xL + xSignDefault * ref.Ax;

    corr = astro.cr3bp.differentialCorrectionPlanarLyapunov( ...
        x0Guess, opts.vy0Guess, mu, opts.maxIter, opts.tol);
end

if ~corr.converged
    error('Planar Lyapunov correction did not converge.');
end

orbit = astro.periodic.buildOrbitStruct(corr.state0, corr.period, mu, ...
    struct('family','planar_lyapunov', ...
           'libr', libr, ...
           'system', 'CR3BP', ...
           'source', ref.source, ...
           'dimension', '2D'), ...
    struct('verbose', opts.verbose));

x = orbit.x(:,1);
y = orbit.x(:,2);

measured = struct();
measured.x0 = corr.state0(1);
measured.vy0 = corr.state0(5);
measured.period = corr.period;
measured.jacobi = astro.cr3bp.jacobiConstant(corr.state0.', mu);
measured.Ax = max(abs(x - xL));
measured.Ay = max(abs(y));
measured.closureError = orbit.closureError;

errors = struct();
errors.Ax = localDiff(ref, measured, 'Ax');
errors.Ay = localDiff(ref, measured, 'Ay');
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
[names, targets, values, absErrs] = localAppend(names, targets, values, absErrs, ref, measured, errors, 'Ay');
[names, targets, values, absErrs] = localAppend(names, targets, values, absErrs, ref, measured, errors, 'period');
[names, targets, values, absErrs] = localAppend(names, targets, values, absErrs, ref, measured, errors, 'jacobi');

if isempty(names)
    out.summaryTable = table();
else
    out.summaryTable = table(string(names(:)), targets(:), values(:), absErrs(:), ...
        'VariableNames', {'Quantity','Target','Measured','AbsError'});
end

if opts.verbose
    fprintf('\nPlanar Lyapunov literature validation\n');
    fprintf('-------------------------------------\n');
    fprintf('Source           : %s\n', ref.source);
    fprintf('Libration point  : L%d\n', libr);
    fprintf('xL               : %.16f\n', xL);
    fprintf('x0 corrected     : %.16f\n', measured.x0);
    fprintf('vy0 corrected    : %.16f\n', measured.vy0);
    fprintf('Period measured  : %.16f\n', measured.period);
    fprintf('Jacobi measured  : %.16f\n', measured.jacobi);
    fprintf('Ax measured      : %.16e\n', measured.Ax);
    fprintf('Ay measured      : %.16e\n', measured.Ay);
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