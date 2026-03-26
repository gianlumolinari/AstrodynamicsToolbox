function eigData = computeEigenstructure(orbit, mu, opts)
%COMPUTEEIGENSTRUCTURE Compute stable/unstable eigendirections along a periodic orbit.
%
% INPUTS
%   orbit : struct with fields .x, .t, .T
%   mu    : CR3BP mass parameter
%   opts  : optional struct
%
% OUTPUT
%   eigData : struct

    if nargin < 3 || isempty(opts)
        opts = struct();
    end

    if ~isfield(opts,'normalizeMode'), opts.normalizeMode = 'fullstate'; end
    if ~isfield(opts,'verbose'),       opts.verbose = true; end
    if ~isfield(opts,'RelTol'),        opts.RelTol = 1e-12; end
    if ~isfield(opts,'AbsTol'),        opts.AbsTol = 1e-12; end

    xOrbit = orbit.x;
    tOrbit = orbit.t;

    if isfield(orbit,'monodromy') && ~isempty(orbit.monodromy)
        M = orbit.monodromy;
    else
        prop0 = astro.cr3bp.propagateWithSTM( ...
            xOrbit(1,:).', [0 orbit.T], mu, ...
            odeset('RelTol',opts.RelTol,'AbsTol',opts.AbsTol));
        M = prop0.PhiFinal;
    end

    [V,D] = eig(M);
    lambda = diag(D);
    mag = abs(lambda);

    if opts.verbose
        fprintf('\nMonodromy eigenvalues:\n');
        for i = 1:numel(lambda)
            fprintf('  lambda(%d) = %+ .12e %+ .12ei   |lambda| = %.12e\n', ...
                i, real(lambda(i)), imag(lambda(i)), abs(lambda(i)));
        end
        fprintf('\n');
    end

    % Select dominant expanding / contracting pair.
    [maxMag, idxUnstable] = max(mag);
    [minMag, idxStable]   = min(mag);

    lambdaUnstable = lambda(idxUnstable);
    lambdaStable   = lambda(idxStable);

    if opts.verbose
        fprintf('Selected unstable eigenvalue: %+ .12e %+ .12ei   |lambda| = %.12e\n', ...
            real(lambdaUnstable), imag(lambdaUnstable), abs(lambdaUnstable));
        fprintf('Selected stable   eigenvalue: %+ .12e %+ .12ei   |lambda| = %.12e\n', ...
            real(lambdaStable), imag(lambdaStable), abs(lambdaStable));
    end

    % This orbit should be hyperbolic if it is a Lyapunov orbit around L1/L2.
    if maxMag < 1 + 1e-6 || minMag > 1 - 1e-6
        error(['Monodromy spectrum is not hyperbolic enough for manifold generation. ' ...
               'This usually means the orbit corrector converged to the wrong closed orbit ' ...
               'or the STM/period is inconsistent.']);
    end

    vUnstable0 = real(V(:,idxUnstable));
    vStable0   = real(V(:,idxStable));

    if norm(vUnstable0) < 1e-14 || norm(vStable0) < 1e-14
        error('Selected eigenvectors are near zero.');
    end

    vUnstable0 = localNormalizeVec(vUnstable0, opts.normalizeMode);
    vStable0   = localNormalizeVec(vStable0, opts.normalizeMode);

    n = size(xOrbit,1);
    vStableAlongOrbit   = zeros(n,6);
    vUnstableAlongOrbit = zeros(n,6);

    x0 = xOrbit(1,:).';

    for k = 1:n
        tk = tOrbit(k);

        if abs(tk) < 1e-14
            Phi = eye(6);
        else
            propSTM = astro.cr3bp.propagateWithSTM( ...
                x0, [0 tk], mu, ...
                odeset('RelTol',opts.RelTol,'AbsTol',opts.AbsTol));

            Phi = propSTM.PhiFinal;
        end

        vs = Phi * vStable0;
        vu = Phi * vUnstable0;

        vs = localNormalizeVec(vs, opts.normalizeMode);
        vu = localNormalizeVec(vu, opts.normalizeMode);

        if k > 1
            if dot(vs, vStableAlongOrbit(k-1,:).') < 0
                vs = -vs;
            end
            if dot(vu, vUnstableAlongOrbit(k-1,:).') < 0
                vu = -vu;
            end
        end

        vStableAlongOrbit(k,:)   = vs.';
        vUnstableAlongOrbit(k,:) = vu.';
    end

    eigData = struct();
    eigData.monodromy           = M;
    eigData.eigvals             = lambda;
    eigData.eigvecs             = V;
    eigData.idxStable           = idxStable;
    eigData.idxUnstable         = idxUnstable;
    eigData.lambdaStable        = lambdaStable;
    eigData.lambdaUnstable      = lambdaUnstable;
    eigData.vStable0            = vStable0;
    eigData.vUnstable0          = vUnstable0;
    eigData.vStableAlongOrbit   = vStableAlongOrbit;
    eigData.vUnstableAlongOrbit = vUnstableAlongOrbit;
end

function v = localNormalizeVec(v, mode)
    switch lower(mode)
        case 'position'
            s = norm(v(1:3));
        case 'fullstate'
            s = norm(v);
        otherwise
            error('Unknown normalization mode.');
    end

    if s < 1e-14
        error('Cannot normalize near-zero eigenvector.');
    end

    v = v / s;
end