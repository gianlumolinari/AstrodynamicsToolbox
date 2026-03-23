function family = continueFamilyPseudoArc(seed1, seed2, nMembers, ds, mu)
%CONTINUEFAMILYPSEUDOARC Pseudo-arclength continuation for planar Lyapunov family.
%
% INPUTS
%   seed1, seed2 : seed orbit structs with fields .state and .period
%   nMembers     : total desired family members (including seeds)
%   ds           : pseudo-arclength step size in reduced parameter space
%   mu           : CR3BP mass parameter
%
% OUTPUT
%   family : struct array with fields
%       .u           = [x0; vy0; thalf]
%       .state0
%       .period
%       .traj
%       .converged
%       .residual
%       .C

if nMembers < 2
    error('nMembers must be at least 2.');
end

family = struct([]);

opts.RelTol = 1e-12;
opts.AbsTol = 1e-12;
opts.Solver = 'ode113';

% ----------------------------------------------------------
% Seed 1
% ----------------------------------------------------------
u1 = [seed1.state(1); seed1.state(5); seed1.period/2];
X01 = [u1(1); 0; 0; 0; u1(2); 0];

traj1 = astro.propagators.propagate( ...
    @(t,x) astro.cr3bp.eomCR3BP(t,x,mu), ...
    [0, 2*u1(3)], X01, opts);

family(1).u = u1;
family(1).state0 = X01;
family(1).period = 2*u1(3);
family(1).traj = traj1.x;
family(1).converged = true;
family(1).residual = 0;
family(1).C = astro.cr3bp.jacobiConstant(X01.', mu);

% ----------------------------------------------------------
% Seed 2
% ----------------------------------------------------------
u2 = [seed2.state(1); seed2.state(5); seed2.period/2];
X02 = [u2(1); 0; 0; 0; u2(2); 0];

traj2 = astro.propagators.propagate( ...
    @(t,x) astro.cr3bp.eomCR3BP(t,x,mu), ...
    [0, 2*u2(3)], X02, opts);

family(2).u = u2;
family(2).state0 = X02;
family(2).period = 2*u2(3);
family(2).traj = traj2.x;
family(2).converged = true;
family(2).residual = 0;
family(2).C = astro.cr3bp.jacobiConstant(X02.', mu);

fprintf('\nStarting pseudo-arclength continuation...\n');
fprintf('Stored seed members 1 and 2\n');

currentDs = ds;
maxRetries = 8;
minDs = max(abs(ds)*1e-4, 1e-6);

for k = 3:nMembers
    fprintf('\nAttempting family member %d\n', k);

    uPrev = family(k-1).u(:);
    uPrev2 = family(k-2).u(:);

    tangent = uPrev - uPrev2;
    tangent = tangent / norm(tangent);

    success = false;
    dsTry = currentDs;

    for retry = 1:maxRetries
        fprintf('  Retry %d with ds = %.3e\n', retry, dsTry);

        uPred = uPrev + dsTry * tangent;

        try
            corr = astro.cr3bp.correctPlanarLyapunovPseudoArc( ...
                uPred, tangent, mu, 20, 1e-11);
        catch ME
            fprintf('    Corrector failed: %s\n', ME.message);
            dsTry = dsTry / 2;
            if abs(dsTry) < minDs
                break
            end
            continue
        end

        if corr.converged
            success = true;
            break
        else
            fprintf('    Corrector did not converge. Residual = %.3e\n', corr.residual);
            dsTry = dsTry / 2;
            if abs(dsTry) < minDs
                break
            end
        end
    end

    if ~success
        fprintf('\nPseudo-arclength continuation stopped at member %d.\n', k);
        fprintf('Minimum ds reached.\n');
        break
    end

    X0 = corr.state0;

    traj = astro.propagators.propagate( ...
        @(t,x) astro.cr3bp.eomCR3BP(t,x,mu), ...
        [0, corr.period], X0, opts);

    family(k).u = corr.u;
    family(k).state0 = X0;
    family(k).period = corr.period;
    family(k).traj = traj.x;
    family(k).converged = true;
    family(k).residual = corr.residual;
    family(k).C = astro.cr3bp.jacobiConstant(X0.', mu);

    fprintf('  Stored family member %d\n', k);
    fprintf('    x0         = %.12f\n', family(k).u(1));
    fprintf('    vy0        = %.12f\n', family(k).u(2));
    fprintf('    halfPeriod = %.12f TU\n', family(k).u(3));
    fprintf('    C          = %.12f\n', family(k).C);

    % Mildly regrow step, but do not exceed original magnitude
    currentDs = sign(ds) * min(abs(ds), 1.2*abs(dsTry));
end

fprintf('\nPseudo-arclength continuation complete. Stored %d members.\n', numel(family));

end