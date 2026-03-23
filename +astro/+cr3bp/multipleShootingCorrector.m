function out = multipleShootingCorrector(Xnodes, dtSeg, mu, maxIter, tol)
%MULTIPLESHOOTINGCORRECTOR Periodic multiple-shooting corrector for CR3BP.
%
% The first node is kept fixed as a phase anchor.
% Nodes 2..N are corrected so that segment continuity and periodic closure
% are restored.
%
% INPUTS
%   Xnodes : [6 x N] node-state guess
%   dtSeg  : scalar or [1 x N] segment durations
%   mu     : CR3BP mass parameter
%   maxIter: maximum Newton iterations
%   tol    : continuity-defect tolerance
%
% OUTPUT
%   out : struct with fields
%       .Xnodes
%       .dtSeg
%       .converged
%       .iterations
%       .residual
%       .segments

if nargin < 4 || isempty(maxIter)
    maxIter = 20;
end
if nargin < 5 || isempty(tol)
    tol = 1e-10;
end

[rows, N] = size(Xnodes);
if rows ~= 6
    error('Xnodes must be 6 x N.');
end

if isscalar(dtSeg)
    dtSeg = repmat(dtSeg, 1, N);
else
    dtSeg = dtSeg(:).';
end

if numel(dtSeg) ~= N
    error('dtSeg must be scalar or length N.');
end

converged = false;
residual = NaN;
segments = cell(N,1);

for k = 1:maxIter
    % Variables are nodes 2..N only
    nVar = 6*(N-1);

    % Constraints are all 6N continuity conditions
    F = zeros(6*N, 1);
    DF = zeros(6*N, nVar);

    for i = 1:N
        Xi = Xnodes(:,i);

        prop = astro.cr3bp.propagateWithSTM( ...
            Xi, [0 dtSeg(i)], mu, odeset('RelTol',1e-12,'AbsTol',1e-12));

        Xf = prop.x(end,:).';
        STM = prop.PhiFinal;

        if i < N
            Xnext = Xnodes(:,i+1);
        else
            Xnext = Xnodes(:,1);
        end

        rowIdx = (6*(i-1)+1):(6*i);
        F(rowIdx) = Xf - Xnext;

        % Derivative wrt current node Xi, if Xi is free (i >= 2)
        if i >= 2
            colIdx_i = (6*(i-2)+1):(6*(i-1));
            DF(rowIdx, colIdx_i) = STM;
        end

        % Derivative wrt next node, if Xnext is free
        % Xnext is free for i = 1,...,N-1 with next >= 2
        if i < N
            colIdx_next = (6*(i-1)+1):(6*i);
            DF(rowIdx, colIdx_next) = DF(rowIdx, colIdx_next) - eye(6);
        else
            % For last segment closure, Xnext = X1 is fixed, so no column
        end

        segments{i}.t = prop.t;
        segments{i}.x = prop.x;
        segments{i}.PhiFinal = STM;
    end

    residual = norm(F);
    fprintf('Multiple shooting iter %2d: ||F|| = %.3e\n', k, residual);

    if residual < tol
        converged = true;
        break
    end

    % Overdetermined least-squares correction
    delta = -DF \ F;

    Xcorr = reshape(delta, 6, N-1);
    Xnodes(:,2:end) = Xnodes(:,2:end) + Xcorr;
end

out.Xnodes = Xnodes;
out.dtSeg = dtSeg;
out.converged = converged;
out.iterations = k;
out.residual = residual;
out.segments = segments;
end