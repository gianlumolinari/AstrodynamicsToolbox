classdef testMultipleShootingCorrectorEarthMoon < matlab.unittest.TestCase
    %TESTMULTIPLESHOOTINGCORRECTOREARTHMOON
    % Validate CR3BP multiple-shooting correction on Earth-Moon halo cases.
    %
    % This is the main robust periodic-correction test for halo orbits.

    methods (Test)
        function testEarthMoonL1HaloNorthMultipleShooting(testCase)
            benchmark = localLoadBenchmark('jpl_em_l1_halo_north.mat');
            mu = localGetMuFromBenchmark(benchmark);
            localRunMultipleShootingTest(testCase, benchmark, mu);
        end

        function testEarthMoonL2HaloNorthMultipleShooting(testCase)
            benchmark = localLoadBenchmark('jpl_em_l2_halo_north.mat');
            mu = localGetMuFromBenchmark(benchmark);
            localRunMultipleShootingTest(testCase, benchmark, mu);
        end
    end
end

function localRunMultipleShootingTest(testCase, benchmark, mu)

    rng(1);

    T = benchmark.orbits;
    n = height(T);

    candidateRows = unique([ ...
        max(2, round(0.30*n)), ...
        max(2, round(0.45*n)), ...
        max(2, round(0.60*n)) ...
    ]);

    success = false;
    bestResidual = inf;
    bestMaxDefect = inf;
    bestRow = NaN;
    bestMessage = '';

    for kk = 1:numel(candidateRows)
        i = candidateRows(kk);

        x0 = [T.x(i); T.y(i); T.z(i); T.vx(i); T.vy(i); T.vz(i)];
        period = T.period(i);

        opts = struct();
        opts.RelTol = 1e-12;
        opts.AbsTol = 1e-12;
        opts.Solver = 'ode113';

        N = 4;
        tNodes = linspace(0, period, N+1);

        fullTraj = astro.propagators.propagate( ...
            @(t,x) astro.cr3bp.eomCR3BP(t,x,mu), ...
            [0 period], x0, opts);

        Xnodes = zeros(6, N);
        for j = 1:N
            tj = tNodes(j);
            Xnodes(:,j) = interp1(fullTraj.t, fullTraj.x, tj).';
        end

        % Keep first node fixed, perturb nodes 2..N only
        Xnodes(:,2:end) = Xnodes(:,2:end) + 1e-8*randn(6, N-1);

        dtSeg = period / N;

        try
            out = astro.cr3bp.multipleShootingCorrector(Xnodes, dtSeg, mu, 20, 1e-10);
        catch ME
            fprintf('\nMultiple shooting candidate row %d failed with error:\n%s\n', i, ME.message);
            bestMessage = ME.message;
            continue
        end

        Xcorr = out.Xnodes;
        maxDefect = 0;

        for j = 1:N
            prop = astro.cr3bp.propagateWithSTM( ...
                Xcorr(:,j), [0 dtSeg], mu, odeset('RelTol',1e-12,'AbsTol',1e-12));

            Xf = prop.x(end,:).';

            if j < N
                defect = norm(Xf - Xcorr(:,j+1));
            else
                defect = norm(Xf - Xcorr(:,1));
            end

            maxDefect = max(maxDefect, defect);
        end

        fprintf('\nMultiple shooting candidate:\n');
        fprintf('  row       = %d\n', i);
        fprintf('  converged = %d\n', out.converged);
        fprintf('  residual  = %.3e\n', out.residual);
        fprintf('  maxDefect = %.3e\n', maxDefect);

        if out.residual < bestResidual
            bestResidual = out.residual;
            bestMaxDefect = maxDefect;
            bestRow = i;
        end

        if out.converged && out.residual < 1e-8 && maxDefect < 1e-8
            success = true;
            break
        end
    end

    fprintf('\nBest multiple-shooting attempt summary:\n');
    fprintf('  best row       = %d\n', bestRow);
    fprintf('  best residual  = %.3e\n', bestResidual);
    fprintf('  best maxDefect = %.3e\n', bestMaxDefect);

    testCase.verifyTrue(success, sprintf([ ...
        'Multiple-shooting corrector did not succeed on any candidate row. ', ...
        'Best row = %d, residual = %.3e, maxDefect = %.3e. %s'], ...
        bestRow, bestResidual, bestMaxDefect, bestMessage));
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