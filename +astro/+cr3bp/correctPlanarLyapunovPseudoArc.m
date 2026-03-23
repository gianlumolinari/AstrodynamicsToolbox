function out = correctPlanarLyapunovPseudoArc(uPred, tangent, mu, maxIter, tol)
%CORRECTPLANARLYAPUNOVPSEUDOARC Pseudo-arclength corrector for planar Lyapunov orbits.
%
% Unknown vector:
%   u = [x0; vy0; thalf]
%
% Initial state form:
%   X0 = [x0; 0; 0; 0; vy0; 0]
%
% Constraints:
%   1) y(thalf)   = 0
%   2) xd(thalf)  = 0
%   3) tangent' * (u - uPred) = 0    (pseudo-arclength hyperplane)
%
% INPUTS
%   uPred   : predictor [3x1]
%   tangent : secant/tangent direction [3x1], assumed normalized
%   mu      : CR3BP mass parameter
%   maxIter : maximum Newton iterations
%   tol     : tolerance on norm(F)
%
% OUTPUT
%   out : struct
%       .u
%       .state0
%       .halfPeriod
%       .period
%       .converged
%       .iterations
%       .residual
%       .trajHalf

if nargin < 4 || isempty(maxIter)
    maxIter = 20;
end
if nargin < 5 || isempty(tol)
    tol = 1e-11;
end

u = uPred(:);
tangent = tangent(:);
tangent = tangent / norm(tangent);

converged = false;
trajHalf = struct('t', [], 'x', []);

for k = 1:maxIter
    x0 = u(1);
    vy0 = u(2);
    thalf = u(3);

    if thalf <= 0
        error('Half-period became nonpositive during pseudo-arclength correction.');
    end

    X0 = [x0; 0; 0; 0; vy0; 0];

    prop = astro.cr3bp.propagateWithSTM( ...
        X0, [0 thalf], mu, odeset('RelTol',1e-12,'AbsTol',1e-12));

    Xf = prop.x(end,:).';
    Phi = prop.PhiFinal;

    % Final dynamics for time derivatives wrt final time
    fFinal = astro.cr3bp.eomCR3BP(0, Xf, mu);

    % Constraints
    F = [ ...
        Xf(2); ...
        Xf(4); ...
        tangent.' * (u - uPred)];

    res = norm(F);

    trajHalf.t = prop.t;
    trajHalf.x = prop.x;

    fprintf('    Pseudo-arc Newton iter %2d: ||F|| = %.3e\n', k, res);

    if res < tol
        converged = true;
        break
    end

    % Jacobian wrt [x0, vy0, thalf]
    J = [ ...
        Phi(2,1), Phi(2,5), fFinal(2); ...
        Phi(4,1), Phi(4,5), fFinal(4); ...
        tangent(1), tangent(2), tangent(3)];

    delta = -J \ F;
    u = u + delta;
end

out.u = u;
out.state0 = [u(1); 0; 0; 0; u(2); 0];
out.halfPeriod = u(3);
out.period = 2*u(3);
out.converged = converged;
out.iterations = k;
out.residual = norm(F);
out.trajHalf = trajHalf;
end