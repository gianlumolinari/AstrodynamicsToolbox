function out = differentialCorrectionHalo(x0, z0, vy0, thalf, mu, maxIter, tol)
%DIFFERENTIALCORRECTIONHALO Differential correction for a 3D halo orbit.
%
% Initial state form:
%   X0 = [x0; 0; z0; 0; vy0; 0]
%
% Unknowns corrected:
%   u = [x0; vy0; thalf]
%
% Fixed parameter:
%   z0
%
% Half-period symmetry constraints at t = thalf:
%   y(thalf)   = 0
%   xd(thalf)  = 0
%   zd(thalf)  = 0
%
% INPUTS
%   x0, z0, vy0, thalf : initial guess parameters
%   mu                 : CR3BP mass parameter
%   maxIter            : maximum Newton iterations
%   tol                : convergence tolerance on norm(F)
%
% OUTPUT
%   out : struct with fields
%       .state0
%       .x0
%       .z0
%       .vy0
%       .halfPeriod
%       .period
%       .converged
%       .iterations
%       .residual
%       .stateHalf
%       .trajHalf
%
% NOTES
%   This is the standard single-shooting symmetry correction:
%       F = [ y_f ; xd_f ; zd_f ] = 0
%   with unknowns [x0, vy0, thalf].

if nargin < 6 || isempty(maxIter)
    maxIter = 20;
end
if nargin < 7 || isempty(tol)
    tol = 1e-11;
end

u = [x0; vy0; thalf];
converged = false;
trajHalf = struct('t', [], 'x', []);

for k = 1:maxIter
    x0 = u(1);
    vy0 = u(2);
    thalf = u(3);

    if thalf <= 0
        error('Half-period became nonpositive during halo correction.');
    end

    X0 = [x0; 0; z0; 0; vy0; 0];

    prop = astro.cr3bp.propagateWithSTM( ...
        X0, [0 thalf], mu, odeset('RelTol',1e-12,'AbsTol',1e-12));

    Xf = prop.x(end,:).';
    Phi = prop.PhiFinal;
    fFinal = astro.cr3bp.eomCR3BP(0, Xf, mu);

    % Constraints
    F = [Xf(2); Xf(4); Xf(6)];
    res = norm(F);

    trajHalf.t = prop.t;
    trajHalf.x = prop.x;

    fprintf('Halo correction iter %2d: ||F|| = %.3e\n', k, res);

    if res < tol
        converged = true;
        break
    end

    % Jacobian wrt [x0, vy0, thalf]
    J = [ ...
        Phi(2,1), Phi(2,5), fFinal(2); ...
        Phi(4,1), Phi(4,5), fFinal(4); ...
        Phi(6,1), Phi(6,5), fFinal(6)];

    delta = -J \ F;
    u = u + delta;
end

X0corr = [u(1); 0; z0; 0; u(2); 0];

out = struct();
out.state0 = X0corr;
out.x0 = u(1);
out.z0 = z0;
out.vy0 = u(2);
out.halfPeriod = u(3);
out.period = 2*u(3);
out.converged = converged;
out.iterations = k;
out.residual = norm(F);
out.stateHalf = Xf;
out.trajHalf = trajHalf;
end