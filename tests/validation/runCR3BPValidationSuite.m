function results = runCR3BPValidationSuite()
%RUNCR3BPVALIDATIONSUITE Run the main CR3BP validation tests.

    thisFile = mfilename('fullpath');
    validationDir = fileparts(thisFile);
    testsDir = fileparts(validationDir);
    repoRoot = fileparts(testsDir);

    cd(repoRoot);
    startup

    testFiles = { ...
        fullfile(repoRoot,'tests','validation','testJplCatalogData.m')
        fullfile(repoRoot,'tests','unit','cr3bp','testJplJacobiConsistency.m')
        fullfile(repoRoot,'tests','unit','cr3bp','testJplJacobiDriftEarthMoon.m')
        fullfile(repoRoot,'tests','unit','cr3bp','testJplHaloClosureEarthMoon.m')
        fullfile(repoRoot,'tests','integration','cr3bp','testContinuationPseudoArcEarthMoon.m')
        fullfile(repoRoot,'tests','unit','cr3bp','testSingleShootingCorrectorEarthMoon.m')
        fullfile(repoRoot,'tests','integration','cr3bp','testMultipleShootingCorrectorEarthMoon.m')
        fullfile(repoRoot,'tests','unit','cr3bp','testMonodromyMatrixEarthMoon.m')
        fullfile(repoRoot,'tests','integration','cr3bp','testLyapunovFamilyTrendEarthMoon.m')
    };

    resultCells = cell(numel(testFiles),1);

    fprintf('\n============================================================\n');
    fprintf('Running CR3BP validation suite\n');
    fprintf('============================================================\n');

    for k = 1:numel(testFiles)
        fprintf('\n------------------------------------------------------------\n');
        fprintf('Running: %s\n', testFiles{k});
        fprintf('------------------------------------------------------------\n');

        r = runtests(testFiles{k});
        resultCells{k} = r(:);
    end

    results = vertcat(resultCells{:});

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