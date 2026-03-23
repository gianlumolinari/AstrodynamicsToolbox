function out = singleShootingCorrector(x0, tf, mu, freeIdx, isTfFree, maxIter, tol)
%SINGLESHOOTINGCORRECTOR Generic single-shooting periodic corrector.
%
% Solves:
%   F = x(tf) - x0 = 0
%
% by correcting selected initial-state components and optionally tf.
% This implementation is intended for refinement from a good seed orbit.
%
% INPUTS
%   x0       : initial state guess [6x1]
%   tf       : final time guess
%   mu       : CR3BP mass parameter
%   freeIdx  : indices of x0 to correct, e.g. [2 3 4 5 6]
%              (leaving x(1) fixed is a common phase anchor)
%   isTfFree : true/false
%   maxIter  : maximum Newton iterations
%   tol      : tolerance on norm(F)
%
% OUTPUT
%   out : struct with fields
%       .x0
%       .tf
%       .converged
%       .iterations
%       .residual
%       .trajectory
%       .PhiFinal

if nargin < 4 || isempty(freeIdx)
    freeIdx = [2 3 4 5 6];
end
if nargin < 5 || isempty(isTfFree)
    isTfFree = true;
end
if nargin < 6 || isempty(maxIter)
    maxIter = 20;
end
if nargin < 7 || isempty(tol)
    tol = 1e-10;
end

x0 = x0(:);
freeIdx = freeIdx(:).';

converged = false;
residual = NaN;
traj = struct('t', [], 'x', []);
PhiFinal = eye(6);

for k = 1:maxIter
    prop = astro.cr3bp.propagateWithSTM( ...
        x0, [0 tf], mu, odeset('RelTol',1e-12,'AbsTol',1e-12));

    xf = prop.x(end,:).';
    PhiFinal = prop.PhiFinal;

    F = xf - x0;
    residual = norm(F);

    traj.t = prop.t;
    traj.x = prop.x;

    fprintf('Single shooting iter %2d: ||F|| = %.3e\n', k, residual);

    if residual < tol
        converged = true;
        break
    end

    % Build correction matrix
    nFree = numel(freeIdx);
    E = zeros(6, nFree);
    for j = 1:nFree
        E(freeIdx(j), j) = 1;
    end

    D = PhiFinal(:, freeIdx) - E;

    if isTfFree
        fFinal = astro.cr3bp.cr3bpRHS(0, xf, mu);
        D = [D, fFinal];
    end

    % Least-squares / Newton step
    delta = -D \ F;

    x0(freeIdx) = x0(freeIdx) + delta(1:nFree);

    if isTfFree
        tf = tf + delta(end);
        if tf <= 0
            error('singleShootingCorrector produced nonpositive final time.');
        end
    end
end

out.x0 = x0;
out.tf = tf;
out.converged = converged;
out.iterations = k;
out.residual = residual;
out.trajectory = traj;
out.PhiFinal = PhiFinal;
end