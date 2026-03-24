function out = correctHaloPseudoArc(uPred, tangent, mu, maxIter, tol)
%CORRECTHALOPSEUDOARC Pseudo-arclength corrector for 3D halo orbits.
%
% Unknown vector:
%   u = [x0; z0; vy0; thalf]
%
% Initial state:
%   X0 = [x0; 0; z0; 0; vy0; 0]
%
% Constraints at half period:
%   y(thalf)  = 0
%   xd(thalf) = 0
%   zd(thalf) = 0
%
% Pseudo-arclength condition:
%   tangent' * (u - uPred) = 0
%
% INPUTS
%   uPred   : predictor [4x1]
%   tangent : tangent/secant direction [4x1]
%   mu      : CR3BP mass parameter
%   maxIter : max Newton iterations
%   tol     : tolerance on norm(F)
%
% OUTPUT
%   out : struct with fields
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
tangent = tangent(:) / norm(tangent);

converged = false;
trajHalf = struct('t', [], 'x', []);

for k = 1:maxIter
    x0    = u(1);
    z0    = u(2);
    vy0   = u(3);
    thalf = u(4);

    if thalf <= 0
        error('Half-period became nonpositive during halo pseudo-arclength correction.');
    end

    X0 = [x0; 0; z0; 0; vy0; 0];

    prop = astro.cr3bp.propagateWithSTM( ...
        X0, [0 thalf], mu, odeset('RelTol',1e-12,'AbsTol',1e-12));

    Xf = prop.x(end,:).';
    Phi = prop.PhiFinal;
    fFinal = astro.cr3bp.eomCR3BP(0, Xf, mu);

    F = [ ...
        Xf(2); ...
        Xf(4); ...
        Xf(6); ...
        tangent.' * (u - uPred)];

    res = norm(F);
    trajHalf.t = prop.t;
    trajHalf.x = prop.x;

    fprintf('    Halo pseudo-arc iter %2d: ||F|| = %.3e\n', k, res);

    if res < tol
        converged = true;
        break
    end

    % Jacobian wrt [x0, z0, vy0, thalf]
    J = [ ...
        Phi(2,1), Phi(2,3), Phi(2,5), fFinal(2); ...
        Phi(4,1), Phi(4,3), Phi(4,5), fFinal(4); ...
        Phi(6,1), Phi(6,3), Phi(6,5), fFinal(6); ...
        tangent(1), tangent(2), tangent(3), tangent(4)];

    delta = -J \ F;
    u = u + delta;
end

out = struct();
out.u = u;
out.state0 = [u(1); 0; u(2); 0; u(3); 0];
out.halfPeriod = u(4);
out.period = 2*u(4);
out.converged = converged;
out.iterations = k;
out.residual = norm(F);
out.trajHalf = trajHalf;
end