classdef testJplJacobiDriftEarthMoon < matlab.unittest.TestCase
    %TESTJPLJACOBIDRIFTEARTHMOON
    % Validate Jacobi conservation while propagating Earth-Moon JPL halo
    % orbits for one catalog period.

    methods (Test)
        function testEarthMoonL1HaloNorthJacobiDrift(testCase)
            benchmark = localLoadBenchmark('jpl_em_l1_halo_north.mat');
            mu = 0.012150585609624;

            localCheckJacobiDrift(testCase, benchmark, mu);
        end

        function testEarthMoonL2HaloNorthJacobiDrift(testCase)
            benchmark = localLoadBenchmark('jpl_em_l2_halo_north.mat');
            mu = 0.012150585609624;

            localCheckJacobiDrift(testCase, benchmark, mu);
        end
    end
end

function localCheckJacobiDrift(testCase, benchmark, mu)

    T = benchmark.orbits;
    vars = lower(string(T.Properties.VariableNames));

    assert(any(vars == "x"),      'Missing x column.');
    assert(any(vars == "y"),      'Missing y column.');
    assert(any(vars == "z"),      'Missing z column.');
    assert(any(vars == "vx"),     'Missing vx column.');
    assert(any(vars == "vy"),     'Missing vy column.');
    assert(any(vars == "vz"),     'Missing vz column.');
    assert(any(vars == "period"), 'Missing period column.');

    n = height(T);
    idx = unique([max(2, round(0.2*n)), max(2, round(0.5*n)), max(2, round(0.8*n))]);

    opts = struct();
    opts.RelTol = 1e-12;
    opts.AbsTol = 1e-12;
    opts.Solver = 'ode113';

    maxDrift = 0;
    worstOrbit = NaN;

    for k = 1:numel(idx)
        i = idx(k);

        x0 = [T.x(i); T.y(i); T.z(i); T.vx(i); T.vy(i); T.vz(i)];
        period = T.period(i);

        out = astro.propagators.propagate( ...
            @(t,x) astro.cr3bp.eomCR3BP(t, x, mu), ...
            [0 period], x0, opts);

        C = astro.cr3bp.jacobiConstant(out.x, mu);
        drift = max(abs(C - C(1)));

        fprintf('row = %d   period = %.15f   Jacobi drift = %.3e\n', ...
            i, period, drift);

        if drift > maxDrift
            maxDrift = drift;
            worstOrbit = i;
        end
    end

    fprintf('Worst Jacobi drift = %.3e at row %d\n', maxDrift, worstOrbit);

    testCase.verifyLessThan(maxDrift, 1e-9, ...
        sprintf('Maximum Jacobi drift too large: %.3e', maxDrift));
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