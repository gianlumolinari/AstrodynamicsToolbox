classdef testJplHaloClosureEarthMoon < matlab.unittest.TestCase
    %TESTJPLHALOCLOSUREEARTHMOON
    % Validate Earth-Moon JPL halo orbits by propagating selected catalog
    % states for one catalog period and checking periodic closure.

    methods (Test)
        function testEarthMoonL1HaloNorthClosure(testCase)
            benchmark = localLoadBenchmark('jpl_em_l1_halo_north.mat');
            mu = 0.012150585609624;

            localCheckClosure(testCase, benchmark, mu);
        end

        function testEarthMoonL2HaloNorthClosure(testCase)
            benchmark = localLoadBenchmark('jpl_em_l2_halo_north.mat');
            mu = 0.012150585609624;

            localCheckClosure(testCase, benchmark, mu);
        end
    end
end

function localCheckClosure(testCase, benchmark, mu)

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

    maxPosErr = 0;
    maxVelErr = 0;
    worstRowPos = NaN;
    worstRowVel = NaN;

    for k = 1:numel(idx)
        i = idx(k);

        x0 = [T.x(i); T.y(i); T.z(i); T.vx(i); T.vy(i); T.vz(i)];
        period = T.period(i);

        out = astro.propagators.propagate( ...
            @(t,x) astro.cr3bp.eomCR3BP(t, x, mu), ...
            [0 period], x0, opts);

        xf = out.x(end,:).';

        posErr = norm(xf(1:3) - x0(1:3));
        velErr = norm(xf(4:6) - x0(4:6));

        fprintf('row = %d   period = %.15f   posErr = %.3e   velErr = %.3e\n', ...
            i, period, posErr, velErr);

        if posErr > maxPosErr
            maxPosErr = posErr;
            worstRowPos = i;
        end

        if velErr > maxVelErr
            maxVelErr = velErr;
            worstRowVel = i;
        end
    end

    fprintf('Worst position closure error = %.3e at row %d\n', maxPosErr, worstRowPos);
    fprintf('Worst velocity closure error = %.3e at row %d\n', maxVelErr, worstRowVel);

    testCase.verifyLessThan(maxPosErr, 1e-8, ...
        sprintf('Maximum position closure error too large: %.3e', maxPosErr));

    testCase.verifyLessThan(maxVelErr, 1e-8, ...
        sprintf('Maximum velocity closure error too large: %.3e', maxVelErr));
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