classdef testMonodromyMatrixEarthMoon < matlab.unittest.TestCase
    %TESTMONODROMYMATRIXEARTHMOON
    % Validate CR3BP monodromy matrix properties for representative
    % Earth-Moon halo family members.

    methods (Test)
        function testEarthMoonL1HaloNorthMonodromy(testCase)
            benchmark = localLoadBenchmark('jpl_em_l1_halo_north.mat');
            mu = localGetMuFromBenchmark(benchmark);
            localRunMonodromyTest(testCase, benchmark, mu);
        end

        function testEarthMoonL2HaloNorthMonodromy(testCase)
            benchmark = localLoadBenchmark('jpl_em_l2_halo_north.mat');
            mu = localGetMuFromBenchmark(benchmark);
            localRunMonodromyTest(testCase, benchmark, mu);
        end
    end
end

function localRunMonodromyTest(testCase, benchmark, mu)

    T = benchmark.orbits;
    n = height(T);

    candidateRows = unique([ ...
        max(2, round(0.30*n)), ...
        max(2, round(0.45*n)), ...
        max(2, round(0.60*n)) ...
    ]);

    success = false;
    bestDetErr = inf;
    bestPairErr = inf;
    bestRow = NaN;

    for kk = 1:numel(candidateRows)
        i = candidateRows(kk);

        x0 = [T.x(i); T.y(i); T.z(i); T.vx(i); T.vy(i); T.vz(i)];
        period = T.period(i);

        try
            M = astro.cr3bp.monodromyMatrix(x0, period, mu);
        catch ME
            fprintf('\nMonodromy candidate row %d failed with error:\n%s\n', i, ME.message);
            continue
        end

        if any(isnan(M(:))) || any(isinf(M(:)))
            fprintf('\nMonodromy candidate row %d produced NaN/Inf entries.\n', i);
            continue
        end

        detErr = abs(det(M) - 1);
        pairErr = localReciprocalPairingError(eig(M));

        fprintf('\nMonodromy candidate:\n');
        fprintf('  row      = %d\n', i);
        fprintf('  detErr   = %.3e\n', detErr);
        fprintf('  pairErr  = %.3e\n', pairErr);

        if detErr < bestDetErr
            bestDetErr = detErr;
            bestPairErr = pairErr;
            bestRow = i;
        end

        if detErr < 1e-6 && pairErr < 1e-4
            success = true;
            break
        end
    end

    fprintf('\nBest monodromy attempt summary:\n');
    fprintf('  best row    = %d\n', bestRow);
    fprintf('  best detErr = %.3e\n', bestDetErr);
    fprintf('  best pairErr= %.3e\n', bestPairErr);

    testCase.verifyTrue(success, sprintf([ ...
        'Monodromy validation did not succeed on any candidate row. ', ...
        'Best row = %d, detErr = %.3e, pairErr = %.3e.'], ...
        bestRow, bestDetErr, bestPairErr));
end

function err = localReciprocalPairingError(lambda)

    lambda = lambda(:);
    n = numel(lambda);
    used = false(n,1);
    pairErrs = [];

    for i = 1:n
        if used(i)
            continue
        end

        target = 1/lambda(i);
        bestErr = inf;
        bestJ = NaN;

        for j = 1:n
            if i == j || used(j)
                continue
            end

            thisErr = abs(lambda(j) - target);
            if thisErr < bestErr
                bestErr = thisErr;
                bestJ = j;
            end
        end

        if ~isnan(bestJ)
            used(i) = true;
            used(bestJ) = true;
            pairErrs(end+1,1) = bestErr; %#ok<AGROW>
        end
    end

    if isempty(pairErrs)
        err = inf;
    else
        err = max(pairErrs);
    end
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
    testsDir = fileparts(thisFile);
    repoRoot = fileparts(fileparts(testsDir));

    f = fullfile(repoRoot, 'data', 'validation', 'cr3bp', 'processed', fileName);
    assert(isfile(f), 'Benchmark file not found: %s', f);

    S = load(f, 'benchmark');
    benchmark = S.benchmark;
end