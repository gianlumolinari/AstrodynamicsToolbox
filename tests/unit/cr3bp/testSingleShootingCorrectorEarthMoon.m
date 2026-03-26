classdef testSingleShootingCorrectorEarthMoon < matlab.unittest.TestCase
    %TESTSINGLESHOOTINGCORRECTOREARTHMOON
    % Validate the generic CR3BP single-shooting corrector on one
    % robust Earth-Moon Lyapunov case only.
    
    methods (Test)
        function testEarthMoonL1LyapunovSingleShooting(testCase)
            benchmark = localLoadBenchmark('jpl_em_l1_lyapunov.mat');
            mu = localGetMuFromBenchmark(benchmark);
            localRunSingleShootingTest(testCase, benchmark, mu);
        end
    end
end

function localRunSingleShootingTest(testCase, benchmark, mu)

    T = benchmark.orbits;
    n = height(T);

    candidateRows = unique([ ...
        max(2, round(0.30*n)), ...
        max(2, round(0.40*n)), ...
        max(2, round(0.50*n)) ...
    ]);

    success = false;
    bestResidual = inf;
    bestMessage = '';
    bestRow = NaN;
    bestPosErr = NaN;
    bestVelErr = NaN;
    bestdC = NaN;

    for kk = 1:numel(candidateRows)
        i = candidateRows(kk);

        xRef = [T.x(i); T.y(i); T.z(i); T.vx(i); T.vy(i); T.vz(i)];
        tfRef = T.period(i);

        xGuess = xRef;
        xGuess(2) = xGuess(2) + 1e-10;
        xGuess(4) = xGuess(4) + 1e-10;
        xGuess(5) = xGuess(5) - 1e-10;
        tfGuess = tfRef*(1 + 1e-10);

        freeIdx = [2 3 4 5 6];
        isTfFree = true;
        maxIter = 25;
        tol = 1e-10;

        try
            out = astro.cr3bp.singleShootingCorrector( ...
                xGuess, tfGuess, mu, freeIdx, isTfFree, maxIter, tol);
        catch ME
            fprintf('\nSingle shooting candidate row %d failed with error:\n%s\n', i, ME.message);
            bestMessage = ME.message;
            continue
        end

        opts = struct();
        opts.RelTol = 1e-12;
        opts.AbsTol = 1e-12;
        opts.Solver = 'ode113';

        prop = astro.propagators.propagate( ...
            @(t,x) astro.cr3bp.eomCR3BP(t,x,mu), ...
            [0 out.tf], out.x0, opts);

        xf = prop.x(end,:).';
        defect = xf - out.x0;

        posErr = norm(defect(1:3));
        velErr = norm(defect(4:6));

        Cref = astro.cr3bp.jacobiConstant(xRef.', mu);
        Ccorr = astro.cr3bp.jacobiConstant(out.x0.', mu);
        dC = abs(Ccorr - Cref);

        fprintf('\nSingle shooting candidate:\n');
        fprintf('  row        = %d\n', i);
        fprintf('  converged  = %d\n', out.converged);
        fprintf('  iterations = %d\n', out.iterations);
        fprintf('  residual   = %.3e\n', out.residual);
        fprintf('  posErr     = %.3e\n', posErr);
        fprintf('  velErr     = %.3e\n', velErr);
        fprintf('  |dC|       = %.3e\n', dC);

        if out.residual < bestResidual
            bestResidual = out.residual;
            bestRow = i;
            bestPosErr = posErr;
            bestVelErr = velErr;
            bestdC = dC;
        end

        if out.converged && out.residual < 1e-8 && posErr < 1e-7 && velErr < 1e-7 && dC < 1e-7
            success = true;
            break
        end
    end

    fprintf('\nBest single-shooting attempt summary:\n');
    fprintf('  best row      = %d\n', bestRow);
    fprintf('  best residual = %.3e\n', bestResidual);
    fprintf('  best posErr   = %.3e\n', bestPosErr);
    fprintf('  best velErr   = %.3e\n', bestVelErr);
    fprintf('  best |dC|     = %.3e\n', bestdC);

    testCase.verifyTrue(success, sprintf([ ...
        'Single-shooting corrector did not succeed on any candidate row. ', ...
        'Best row = %d, residual = %.3e, posErr = %.3e, velErr = %.3e, |dC| = %.3e. %s'], ...
        bestRow, bestResidual, bestPosErr, bestVelErr, bestdC, bestMessage));
end

function mu = localGetMuFromBenchmark(benchmark)
    if isfield(benchmark,'raw') && isfield(benchmark.raw,'system') && isfield(benchmark.raw.system,'mass_ratio')
        mu = str2double(string(benchmark.raw.system.mass_ratio));
    elseif isfield(benchmark,'systemInfo') && isfield(benchmark.systemInfo,'mass_ratio')
        mu = str2double(string(benchmark.systemInfo.mass_ratio));
    else
        error('Could not extract mass_ratio from benchmark metadata.');
    end
end

function benchmark = localLoadBenchmark(fileName)
    thisFile = mfilename('fullpath');
    testDir = fileparts(thisFile);
    repoRoot = fileparts(fileparts(fileparts(testDir)));

    f = fullfile(repoRoot, 'data', 'validation', 'cr3bp', 'processed', fileName);
    assert(isfile(f), 'Benchmark file not found: %s', f);

    S = load(f, 'benchmark');
    benchmark = S.benchmark;
end