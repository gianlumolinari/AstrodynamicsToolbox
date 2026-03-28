function report = validatePeriodicOrbitForManifolds(orbit, mu, opts)
%VALIDATEPERIODICORBITFORMANIFOLDS Basic validation of a periodic orbit
% before manifold generation.
%
% INPUTS
%   orbit : struct with fields
%       .state0
%       .T
%       .x
%       .t
%       .monodromy   (optional)
%       .closureError (optional)
%   mu    : CR3BP mass parameter
%   opts  : optional struct
%
% OUTPUT
%   report : struct with fields
%       .closureError
%       .jacobiDrift
%       .eigvals
%       .maxEigMod
%       .minEigMod
%       .idxUnstable
%       .idxStable
%       .lambdaUnstable
%       .lambdaStable
%       .hasHyperbolicPair
%       .isAcceptable
%       .messages

if nargin < 3 || isempty(opts)
    opts = struct();
end

if ~isfield(opts,'closureTol'), opts.closureTol = 1e-8; end
if ~isfield(opts,'jacobiTol'),  opts.jacobiTol  = 1e-10; end
if ~isfield(opts,'verbose'),    opts.verbose    = true; end

messages = {};

% -------------------------------------------------------------------------
% Closure error
% -------------------------------------------------------------------------
if isfield(orbit,'closureError') && ~isempty(orbit.closureError)
    closureError = orbit.closureError;
else
    closureError = norm(orbit.x(end,:).' - orbit.state0);
end

if closureError > opts.closureTol
    messages{end+1} = sprintf('Closure error %.3e exceeds tolerance %.3e.', ...
        closureError, opts.closureTol);
end

% -------------------------------------------------------------------------
% Jacobi drift along the periodic orbit
% -------------------------------------------------------------------------
C = zeros(size(orbit.t));
for k = 1:numel(orbit.t)
    C(k) = astro.cr3bp.jacobiConstant(orbit.x(k,:).', mu);
end
jacobiDrift = max(abs(C - C(1)));

if jacobiDrift > opts.jacobiTol
    messages{end+1} = sprintf('Jacobi drift %.3e exceeds tolerance %.3e.', ...
        jacobiDrift, opts.jacobiTol);
end

% -------------------------------------------------------------------------
% Monodromy / eigenvalues
% -------------------------------------------------------------------------
if isfield(orbit,'monodromy') && ~isempty(orbit.monodromy)
    M = orbit.monodromy;
else
    prop = astro.cr3bp.propagateWithSTM( ...
        orbit.state0, [0 orbit.T], mu, odeset('RelTol',1e-12,'AbsTol',1e-12));
    M = prop.PhiFinal;
end

[V,D] = eig(M);
lambda = diag(D);
lamMod = abs(lambda);

[maxEigMod, idxUnstable] = max(lamMod);
[minEigMod, idxStable]   = min(lamMod);

lambdaUnstable = lambda(idxUnstable);
lambdaStable   = lambda(idxStable);

hasHyperbolicPair = (maxEigMod > 1 + 1e-6) && (minEigMod < 1 - 1e-6);

if ~hasHyperbolicPair
    messages{end+1} = 'No clear hyperbolic stable/unstable pair detected.';
end

isAcceptable = isempty(messages);

report = struct();
report.closureError = closureError;
report.jacobiDrift = jacobiDrift;
report.eigvals = lambda;
report.maxEigMod = maxEigMod;
report.minEigMod = minEigMod;
report.idxUnstable = idxUnstable;
report.idxStable = idxStable;
report.lambdaUnstable = lambdaUnstable;
report.lambdaStable = lambdaStable;
report.hasHyperbolicPair = hasHyperbolicPair;
report.isAcceptable = isAcceptable;
report.messages = messages;

if opts.verbose
    fprintf('\nPeriodic-orbit validation report\n');
    fprintf('--------------------------------\n');
    fprintf('Closure error     : %.3e\n', closureError);
    fprintf('Jacobi drift      : %.3e\n', jacobiDrift);
    fprintf('Largest |lambda|  : %.6e\n', maxEigMod);
    fprintf('Smallest |lambda| : %.6e\n', minEigMod);

    if hasHyperbolicPair
        fprintf('Hyperbolic pair   : YES\n');
    else
        fprintf('Hyperbolic pair   : NO\n');
    end

    if isempty(messages)
        fprintf('Status            : ACCEPTABLE FOR MANIFOLDS\n');
    else
        fprintf('Status            : WARNING\n');
        for i = 1:numel(messages)
            fprintf('  - %s\n', messages{i});
        end
    end
end
end