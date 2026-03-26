function results = runCR3BPValidationSuite()
%RUNCR3BPVALIDATIONSUITE Run the main CR3BP validation tests.
%
% OUTPUT
%   results : array of matlab.unittest.TestResult
%
% Usage:
%   results = runCR3BPValidationSuite;

    startup

    testFiles = { ...
        'tests/validation/testJplCatalogData.m'
        'tests/validation/testJplJacobiConsistency.m'
        'tests/validation/testJplJacobiDriftEarthMoon.m'
        'tests/validation/testJplHaloClosureEarthMoon.m'
        'tests/validation/testContinuationPseudoArcEarthMoon.m'
        'tests/validation/testSingleShootingCorrectorEarthMoon.m'
        'tests/validation/testMultipleShootingCorrectorEarthMoon.m'
        'tests/validation/testMonodromyMatrixEarthMoon.m'
        'tests/validation/testLyapunovFamilyTrendEarthMoon.m'
    };

    allResults = matlab.unittest.TestResult.empty(0,1);

    fprintf('\n============================================================\n');
    fprintf('Running CR3BP validation suite\n');
    fprintf('============================================================\n');

    for k = 1:numel(testFiles)
        fprintf('\n------------------------------------------------------------\n');
        fprintf('Running: %s\n', testFiles{k});
        fprintf('------------------------------------------------------------\n');

        r = runtests(testFiles{k});
        r = r(:);   % force column vector

        allResults = [allResults; r]; %#ok<AGROW>
    end

    results = allResults;

    fprintf('\n============================================================\n');
    fprintf('CR3BP validation suite complete\n');
    fprintf('============================================================\n');

    disp(table(results))

    nTotal = numel(results);
    nPassed = sum([results.Passed]);
    nFailed = sum([results.Failed]);
    nIncomplete = sum([results.Incomplete]);

    fprintf('\nSummary:\n');
    fprintf('  Total tests  : %d\n', nTotal);
    fprintf('  Passed       : %d\n', nPassed);
    fprintf('  Failed       : %d\n', nFailed);
    fprintf('  Incomplete   : %d\n', nIncomplete);

    if nFailed == 0 && nIncomplete == 0
        fprintf('\nAll CR3BP validation tests passed.\n');
    else
        fprintf('\nSome CR3BP validation tests did not pass.\n');
    end
end