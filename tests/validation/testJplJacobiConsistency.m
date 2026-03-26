classdef testJplJacobiConsistency < matlab.unittest.TestCase
    %TESTJPLJACOBICONSISTENCY
    % Check that the Jacobi constant recomputed from frozen JPL states
    % matches the Jacobi value stored in the JPL catalog.

    methods (Test)
        function testEarthMoonL1HaloNorth(testCase)
            benchmark = localLoadBenchmark('jpl_em_l1_halo_north.mat');
            mu = localGetMu(benchmark.system);
            localCheckJacobiConsistency(testCase, benchmark, mu);
        end

        function testEarthMoonL2HaloNorth(testCase)
            benchmark = localLoadBenchmark('jpl_em_l2_halo_north.mat');
            mu = localGetMu(benchmark.system);
            localCheckJacobiConsistency(testCase, benchmark, mu);
        end

        function testSunEarthL1HaloNorth(testCase)
            benchmark = localLoadBenchmark('jpl_se_l1_halo_north.mat');
            mu = localGetMu(benchmark.system);
            localCheckJacobiConsistency(testCase, benchmark, mu);
        end

        function testSunEarthL2HaloNorth(testCase)
            benchmark = localLoadBenchmark('jpl_se_l2_halo_north.mat');
            mu = localGetMu(benchmark.system);
            localCheckJacobiConsistency(testCase, benchmark, mu);
        end
    end
end

function localCheckJacobiConsistency(testCase, benchmark, mu)

    T = benchmark.orbits;
    vars = lower(string(T.Properties.VariableNames));

    assert(any(vars == "x"),      'Missing x column.');
    assert(any(vars == "y"),      'Missing y column.');
    assert(any(vars == "z"),      'Missing z column.');
    assert(any(vars == "vx"),     'Missing vx column.');
    assert(any(vars == "vy"),     'Missing vy column.');
    assert(any(vars == "vz"),     'Missing vz column.');
    assert(any(vars == "jacobi"), 'Missing jacobi column.');

    n = height(T);
    idx = unique(round(linspace(1, n, min(25, n))));

    maxErr = 0;

    for k = 1:numel(idx)
        i = idx(k);

        X = [T.x(i), T.y(i), T.z(i), T.vx(i), T.vy(i), T.vz(i)];
        CjCatalog = T.jacobi(i);

        CjComputed = astro.cr3bp.jacobiConstant(X, mu);

        err = abs(CjComputed - CjCatalog);
        maxErr = max(maxErr, err);
    end

    testCase.verifyLessThan(maxErr, 1e-10, ...
        sprintf('Max Jacobi mismatch too large: %.3e', maxErr));
end

function mu = localGetMu(systemName)

    switch lower(string(systemName))
        case "earth-moon"
            mu = 0.012150585609624;

        case "sun-earth"
            mu = 3.054200000000000E-6;

        otherwise
            error('Unsupported CR3BP system: %s', systemName);
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