classdef testContinuationPseudoArcEarthMoon < matlab.unittest.TestCase
    %TESTCONTINUATIONPSEUDOARCEARTHMOON
    % Validate planar Lyapunov pseudo-arclength continuation in the
    % Earth-Moon CR3BP.

    methods (Test)
        function testEarthMoonL1LyapunovPseudoArc(testCase)
            benchmark = localLoadBenchmark('jpl_em_l1_lyapunov.mat');
            mu = 0.012150585609624;

            localRunContinuationTest(testCase, benchmark, mu);
        end

        function testEarthMoonL2LyapunovPseudoArc(testCase)
            benchmark = localLoadBenchmark('jpl_em_l2_lyapunov.mat');
            mu = 0.012150585609624;

            localRunContinuationTest(testCase, benchmark, mu);
        end
    end
end

function localRunContinuationTest(testCase, benchmark, mu)

    T = benchmark.orbits;
    n = height(T);

    i1 = max(2, round(0.35*n));
    i2 = max(i1+1, round(0.40*n));

    seed1.state = [T.x(i1); T.y(i1); T.z(i1); T.vx(i1); T.vy(i1); T.vz(i1)];
    seed1.period = T.period(i1);

    seed2.state = [T.x(i2); T.y(i2); T.z(i2); T.vx(i2); T.vy(i2); T.vz(i2)];
    seed2.period = T.period(i2);

    nMembers = 6;
    ds = 1e-3;

    family = astro.cr3bp.continueFamilyPseudoArc(seed1, seed2, nMembers, ds, mu);

    nBuilt = numel(family);
    testCase.verifyGreaterThanOrEqual(nBuilt, 3, ...
        sprintf('Continuation produced too few members: %d', nBuilt));

    Cvals = nan(nBuilt,1);
    Tvals = nan(nBuilt,1);

    for k = 1:nBuilt
        testCase.verifyTrue(family(k).converged, ...
            sprintf('Family member %d is not marked converged.', k));

        testCase.verifyFalse(any(isnan(family(k).state0)), ...
            sprintf('Family member %d contains NaNs.', k));

        Cvals(k) = family(k).C;
        Tvals(k) = family(k).period;
    end

    dC = diff(Cvals);
    dT = diff(Tvals);

    fprintf('Continuation:\n');
    fprintf('  number of members = %d\n', nBuilt);
    fprintf('  C values =\n');
    disp(Cvals.')
    fprintf('  periods =\n');
    disp(Tvals.')

    testCase.verifyTrue(all(abs(dC) < 1), 'Unphysical jump in Jacobi constant.');
    testCase.verifyTrue(all(abs(dT) < 10), 'Unphysical jump in period.');

    % Check neighboring seed parameter vectors vary smoothly
    for k = 2:nBuilt
        du = norm(family(k).u(:) - family(k-1).u(:));
        testCase.verifyLessThan(du, 1, ...
            sprintf('Continuation step between members %d and %d too large: %.3e', ...
            k-1, k, du));
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