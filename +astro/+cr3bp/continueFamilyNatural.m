function family = continueFamilyNatural(seedOrbit, nMembers, stepSize, nNodes, mu)
%CONTINUEFAMILYNATURAL Natural parameter continuation with adaptive step control.
%
% INPUTS
%   seedOrbit : struct with fields .state and .period
%   nMembers  : desired total number of family members, including seed
%   stepSize  : initial continuation step in x0(1)
%   nNodes    : number of multiple-shooting nodes
%   mu        : CR3BP mass parameter
%
% OUTPUT
%   family : struct array with fields
%       .x0
%       .T
%       .traj
%       .converged
%       .stepUsed
%       .residual
%
% NOTES
%   - Member 1 is the seed orbit
%   - Uses x0(1) as the natural continuation parameter
%   - Retries failed steps with reduced step size
%   - Stores only converged family members

family = struct([]);

opts.RelTol = 1e-12;
opts.AbsTol = 1e-12;
opts.Solver = 'ode113';

maxRetries = 8;
minStep = max(abs(stepSize)*1e-4, 1e-6);

% ----------------------------------------------------------
% Store seed orbit as first family member
% ----------------------------------------------------------
xRef = seedOrbit.state(:);
Tref = seedOrbit.period;

seedOut = astro.propagators.propagate( ...
    @(t,x) astro.cr3bp.eomCR3BP(t,x,mu), ...
    linspace(0, Tref, 800), xRef, opts);

family(1).x0 = xRef;
family(1).T = Tref;
family(1).traj = seedOut.x;
family(1).converged = true;
family(1).stepUsed = 0.0;
family(1).residual = 0.0;

fprintf('\nStarting natural continuation...\n');
fprintf('Stored seed orbit as family member 1\n');
fprintf('  x0(1)    = %.12f\n', xRef(1));
fprintf('  Period   = %.12f TU\n', Tref);

nConv = 1;
currentStep = stepSize;

while nConv < nMembers
    fprintf('\nAttempting family member %d\n', nConv + 1);

    success = false;
    stepTry = currentStep;

    for retry = 1:maxRetries
        fprintf('  Retry %d with step = %.3e\n', retry, stepTry);

        % --------------------------------------------------
        % Predictor: perturb x0(1)
        % --------------------------------------------------
        xPred = xRef;
        xPred(1) = xPred(1) + stepTry;

        % --------------------------------------------------
        % Build node guess by propagating predictor
        % --------------------------------------------------
        tGrid = linspace(0, Tref, nNodes+1);

        try
            predOut = astro.propagators.propagate( ...
                @(t,x) astro.cr3bp.eomCR3BP(t,x,mu), ...
                tGrid, xPred, opts);
        catch ME
            fprintf('    Predictor propagation failed: %s\n', ME.message);
            stepTry = stepTry / 2;
            if abs(stepTry) < minStep
                break
            end
            continue
        end

        Xnodes = predOut.x(1:nNodes,:).';
        dtSeg = (Tref/nNodes) * ones(1,nNodes);

        % --------------------------------------------------
        % Correct with multiple shooting
        % --------------------------------------------------
        try
            out = astro.cr3bp.multipleShootingCorrector( ...
                Xnodes, dtSeg, mu, 20, 1e-10);
        catch ME
            fprintf('    Corrector failed: %s\n', ME.message);
            stepTry = stepTry / 2;
            if abs(stepTry) < minStep
                break
            end
            continue
        end

        if out.converged
            success = true;
            break
        else
            fprintf('    Corrector did not converge. Residual = %.3e\n', out.residual);
            stepTry = stepTry / 2;
            if abs(stepTry) < minStep
                break
            end
        end
    end

    if ~success
        fprintf('\nContinuation stopped: could not converge next family member.\n');
        fprintf('Minimum attempted step size reached.\n');
        break
    end

    % ------------------------------------------------------
    % Reconstruct full trajectory
    % ------------------------------------------------------
    traj = [];
    for i = 1:numel(out.segments)
        Xi = out.segments{i}.x;
        if i < numel(out.segments)
            traj = [traj; Xi(1:end-1,:)]; %#ok<AGROW>
        else
            traj = [traj; Xi]; %#ok<AGROW>
        end
    end

    nConv = nConv + 1;
    family(nConv).x0 = out.Xnodes(:,1);
    family(nConv).T = sum(out.dtSeg);
    family(nConv).traj = traj;
    family(nConv).converged = out.converged;
    family(nConv).stepUsed = stepTry;
    family(nConv).residual = out.residual;

    fprintf('  Stored family member %d\n', nConv);
    fprintf('    x0(1)    = %.12f\n', family(nConv).x0(1));
    fprintf('    Period   = %.12f TU\n', family(nConv).T);
    fprintf('    Residual = %.3e\n', family(nConv).residual);

    % Update reference for next step
    xRef = family(nConv).x0;
    Tref = family(nConv).T;

    % Mildly regrow step after success, but do not exceed original magnitude
    currentStep = sign(stepSize) * min(abs(stepSize), 1.2*abs(stepTry));
end

fprintf('\nContinuation complete. Stored %d converged family members.\n', numel(family));

end