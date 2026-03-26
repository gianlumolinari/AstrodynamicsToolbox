classdef testLyapunovFamilyTrendEarthMoon < matlab.unittest.TestCase
    %TESTLYAPUNOVFAMILYTRENDEARTHMOON
    % Compare pseudo-arclength continuation output against the frozen JPL
    % Lyapunov family trend in a local sense.

    methods (Test)
        function testEarthMoonL1LyapunovTrend(testCase)
            benchmark = localLoadBenchmark('jpl_em_l1_lyapunov.mat');
            mu = localGetMuFromBenchmark(benchmark);
            localRunTrendTest(testCase, benchmark, mu);
        end

        function testEarthMoonL2LyapunovTrend(testCase)
            benchmark = localLoadBenchmark('jpl_em_l2_lyapunov.mat');
            mu = localGetMuFromBenchmark(benchmark);
            localRunTrendTest(testCase, benchmark, mu);
        end
    end
end

function localRunTrendTest(testCase, benchmark, mu)

    T = benchmark.orbits;
    n = height(T);

    i1 = max(2, round(0.35*n));
    i2 = max(i1+1, round(0.40*n));

    seed1.state = [T.x(i1); T.y(i1); T.z(i1); T.vx(i1); T.vy(i1); T.vz(i1)];
    seed1.period = T.period(i1);

    seed2.state = [T.x(i2); T.y(i2); T.z(i2); T.vx(i2); T.vy(i2); T.vz(i2)];
    seed2.period = T.period(i2);

    family = astro.cr3bp.continueFamilyPseudoArc(seed1, seed2, 6, 1e-3, mu);
    nBuilt = numel(family);

    testCase.verifyGreaterThanOrEqual(nBuilt, 3, ...
        sprintf('Continuation produced too few members: %d', nBuilt));

    Cfam = nan(nBuilt,1);
    Pfam = nan(nBuilt,1);
    x0fam = nan(nBuilt,1);
    vy0fam = nan(nBuilt,1);

    for k = 1:nBuilt
        testCase.verifyTrue(family(k).converged, ...
            sprintf('Family member %d is not converged.', k));

        Cfam(k) = family(k).C;
        Pfam(k) = family(k).period;
        x0fam(k) = family(k).state0(1);
        vy0fam(k) = family(k).state0(5);
    end

    % Local JPL reference window around seeds
    jIdx = unique(round(linspace(i1, min(n,i1+40), nBuilt)));
    Cj = T.jacobi(jIdx);
    Pj = T.period(jIdx);
    x0j = T.x(jIdx);
    vy0j = T.vy(jIdx);

    % Compare first differences rather than absolute orbit-by-orbit equality
    dCf = diff(Cfam);
    dPf = diff(Pfam);
    dxF = diff(x0fam);
    dvyF = diff(vy0fam);

    dCj = diff(Cj);
    dPj = diff(Pj);
    dxJ = diff(x0j);
    dvyJ = diff(vy0j);

    fprintf('Lyapunov family trend:\n');
    fprintf('  nBuilt = %d\n', nBuilt);

    testCase.verifyTrue(all(isfinite(Cfam)));
    testCase.verifyTrue(all(isfinite(Pfam)));
    testCase.verifyTrue(all(isfinite(x0fam)));
    testCase.verifyTrue(all(isfinite(vy0fam)));

    testCase.verifyTrue(all(abs(dCf) < 1), 'Unphysical jump in continuation Jacobi.');
    testCase.verifyTrue(all(abs(dPf) < 10), 'Unphysical jump in continuation period.');

    % Trend sign consistency
    testCase.verifyTrue(all(signNonzero(dPf) == signNonzero(dPj(1:numel(dPf)))), ...
        'Period trend sign is inconsistent with JPL.');
    testCase.verifyTrue(all(signNonzero(dxF) == signNonzero(dxJ(1:numel(dxF)))), ...
        'x0 trend sign is inconsistent with JPL.');
    testCase.verifyTrue(all(signNonzero(dvyF) == signNonzero(dvyJ(1:numel(dvyF)))), ...
        'vy0 trend sign is inconsistent with JPL.');

    % Loose magnitude sanity
    testCase.verifyLessThan(max(abs(dPf - dPj(1:numel(dPf)))), 1, ...
        'Continuation period increments deviate too much from JPL.');
end

function s = signNonzero(v)
    s = sign(v);
    s(s == 0) = 1;
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