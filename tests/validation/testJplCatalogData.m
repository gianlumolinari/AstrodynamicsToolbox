classdef testJplCatalogData < matlab.unittest.TestCase
    %TESTJPLCATALOGDATA
    % Smoke tests for frozen JPL CR3BP benchmark files.
    %
    % These tests do NOT propagate trajectories yet.
    % They only verify that:
    %   1) the processed MAT benchmark files exist
    %   2) the saved benchmark struct contains the expected metadata
    %   3) the orbit table has the essential fields needed for later tests

    methods (Test)
        function testEarthMoonL1HaloNorth(testCase)
            benchmark = localLoadBenchmark('jpl_em_l1_halo_north.mat');

            testCase.verifyEqual(benchmark.system, 'earth-moon');
            testCase.verifyEqual(benchmark.family, 'halo');
            testCase.verifyEqual(benchmark.librationPoint, 1);
            testCase.verifyEqual(benchmark.branch, 'N');
            testCase.verifyGreaterThan(benchmark.nOrbits, 0);

            vars = benchmark.orbits.Properties.VariableNames;
            testCase.verifyTrue(any(strcmpi(vars, 'x')),      'Missing field x');
            testCase.verifyTrue(any(strcmpi(vars, 'y')),      'Missing field y');
            testCase.verifyTrue(any(strcmpi(vars, 'z')),      'Missing field z');
            testCase.verifyTrue(any(strcmpi(vars, 'vx')),     'Missing field vx');
            testCase.verifyTrue(any(strcmpi(vars, 'vy')),     'Missing field vy');
            testCase.verifyTrue(any(strcmpi(vars, 'vz')),     'Missing field vz');
            testCase.verifyTrue(any(strcmpi(vars, 'period')), 'Missing field period');
            testCase.verifyTrue(any(strcmpi(vars, 'jacobi')), 'Missing field jacobi');
        end

        function testEarthMoonL2HaloNorth(testCase)
            benchmark = localLoadBenchmark('jpl_em_l2_halo_north.mat');

            testCase.verifyEqual(benchmark.system, 'earth-moon');
            testCase.verifyEqual(benchmark.family, 'halo');
            testCase.verifyEqual(benchmark.librationPoint, 2);
            testCase.verifyEqual(benchmark.branch, 'N');
            testCase.verifyGreaterThan(benchmark.nOrbits, 0);

            vars = benchmark.orbits.Properties.VariableNames;
            testCase.verifyTrue(any(strcmpi(vars, 'x')),      'Missing field x');
            testCase.verifyTrue(any(strcmpi(vars, 'y')),      'Missing field y');
            testCase.verifyTrue(any(strcmpi(vars, 'z')),      'Missing field z');
            testCase.verifyTrue(any(strcmpi(vars, 'vx')),     'Missing field vx');
            testCase.verifyTrue(any(strcmpi(vars, 'vy')),     'Missing field vy');
            testCase.verifyTrue(any(strcmpi(vars, 'vz')),     'Missing field vz');
            testCase.verifyTrue(any(strcmpi(vars, 'period')), 'Missing field period');
            testCase.verifyTrue(any(strcmpi(vars, 'jacobi')), 'Missing field jacobi');
        end

        function testSunEarthL1HaloNorth(testCase)
            benchmark = localLoadBenchmark('jpl_se_l1_halo_north.mat');

            testCase.verifyEqual(benchmark.system, 'sun-earth');
            testCase.verifyEqual(benchmark.family, 'halo');
            testCase.verifyEqual(benchmark.librationPoint, 1);
            testCase.verifyEqual(benchmark.branch, 'N');
            testCase.verifyGreaterThan(benchmark.nOrbits, 0);

            vars = benchmark.orbits.Properties.VariableNames;
            testCase.verifyTrue(any(strcmpi(vars, 'x')),      'Missing field x');
            testCase.verifyTrue(any(strcmpi(vars, 'y')),      'Missing field y');
            testCase.verifyTrue(any(strcmpi(vars, 'z')),      'Missing field z');
            testCase.verifyTrue(any(strcmpi(vars, 'vx')),     'Missing field vx');
            testCase.verifyTrue(any(strcmpi(vars, 'vy')),     'Missing field vy');
            testCase.verifyTrue(any(strcmpi(vars, 'vz')),     'Missing field vz');
            testCase.verifyTrue(any(strcmpi(vars, 'period')), 'Missing field period');
            testCase.verifyTrue(any(strcmpi(vars, 'jacobi')), 'Missing field jacobi');
        end

        function testSunEarthL2HaloNorth(testCase)
            benchmark = localLoadBenchmark('jpl_se_l2_halo_north.mat');

            testCase.verifyEqual(benchmark.system, 'sun-earth');
            testCase.verifyEqual(benchmark.family, 'halo');
            testCase.verifyEqual(benchmark.librationPoint, 2);
            testCase.verifyEqual(benchmark.branch, 'N');
            testCase.verifyGreaterThan(benchmark.nOrbits, 0);

            vars = benchmark.orbits.Properties.VariableNames;
            testCase.verifyTrue(any(strcmpi(vars, 'x')),      'Missing field x');
            testCase.verifyTrue(any(strcmpi(vars, 'y')),      'Missing field y');
            testCase.verifyTrue(any(strcmpi(vars, 'z')),      'Missing field z');
            testCase.verifyTrue(any(strcmpi(vars, 'vx')),     'Missing field vx');
            testCase.verifyTrue(any(strcmpi(vars, 'vy')),     'Missing field vy');
            testCase.verifyTrue(any(strcmpi(vars, 'vz')),     'Missing field vz');
            testCase.verifyTrue(any(strcmpi(vars, 'period')), 'Missing field period');
            testCase.verifyTrue(any(strcmpi(vars, 'jacobi')), 'Missing field jacobi');
        end
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