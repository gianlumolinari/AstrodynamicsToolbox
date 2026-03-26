function buildJplCr3bpValidationData()
%BUILDJPLCR3BPVALIDATIONDATA
% Freeze selected JPL periodic-orbit catalog families into local JSON and MAT
% files for offline validation tests.
%
% This function:
%   1) queries the JPL Three-Body Periodic Orbit API through the existing
%      repo function astro.cr3bp.queryJPLPeriodicOrbits
%   2) saves the raw API response as JSON
%   3) converts the returned family into a clean MATLAB benchmark struct
%   4) saves the benchmark as a .mat file
%
% Output folders:
%   data/validation/cr3bp/raw
%   data/validation/cr3bp/processed
%
% Usage:
%   buildJplCr3bpValidationData

    thisFile = mfilename('fullpath');
    examplesDir = fileparts(thisFile);
    repoRoot = fileparts(fileparts(examplesDir));

    rawDir = fullfile(repoRoot, 'data', 'validation', 'cr3bp', 'raw');
    procDir = fullfile(repoRoot, 'data', 'validation', 'cr3bp', 'processed');

    if ~exist(rawDir, 'dir')
        mkdir(rawDir);
    end
    if ~exist(procDir, 'dir')
        mkdir(procDir);
    end

    cases = [ ...
    struct('sys','earth-moon','family','halo',     'libr',1,'branch','N','tag','jpl_em_l1_halo_north'), ...
    struct('sys','earth-moon','family','halo',     'libr',2,'branch','N','tag','jpl_em_l2_halo_north'), ...
    struct('sys','sun-earth',  'family','halo',    'libr',1,'branch','N','tag','jpl_se_l1_halo_north'), ...
    struct('sys','sun-earth',  'family','halo',    'libr',2,'branch','N','tag','jpl_se_l2_halo_north'), ...
    struct('sys','earth-moon','family','lyapunov', 'libr',1,'branch','','tag','jpl_em_l1_lyapunov'), ...
    struct('sys','earth-moon','family','lyapunov', 'libr',2,'branch','','tag','jpl_em_l2_lyapunov') ...
];
    for k = 1:numel(cases)
        c = cases(k);

        fprintf('\n============================================================\n');
        fprintf('Querying JPL family: %s\n', c.tag);
        fprintf('  sys    = %s\n', c.sys);
        fprintf('  family = %s\n', c.family);
        fprintf('  libr   = %d\n', c.libr);
        fprintf('  branch = %s\n', c.branch);

        data = astro.cr3bp.queryJPLPeriodicOrbits( ...
            'sys', c.sys, ...
            'family', c.family, ...
            'libr', c.libr, ...
            'branch', c.branch);

        rawJsonFile = fullfile(rawDir, [c.tag '.json']);
        procMatFile = fullfile(procDir, [c.tag '.mat']);

        % -----------------------------
        % Save raw JSON response
        % -----------------------------
        jsonText = jsonencode(data);

        fid = fopen(rawJsonFile, 'w');
        if fid == -1
            error('Could not open file for writing: %s', rawJsonFile);
        end
        fwrite(fid, jsonText, 'char');
        fclose(fid);

        % -----------------------------
        % Build processed benchmark
        % -----------------------------
        benchmark = struct();
        benchmark.source = 'JPL Three-Body Periodic Orbits API';
        benchmark.tag = c.tag;
        benchmark.system = c.sys;
        benchmark.family = c.family;
        benchmark.librationPoint = c.libr;
        benchmark.branch = c.branch;
        benchmark.generatedOn = datestr(now, 30);

        % Keep the full response too, for debugging
        benchmark.raw = data;

        % Optional metadata if present in the response
        if isfield(data, 'signature')
            benchmark.signature = data.signature;
        else
            benchmark.signature = [];
        end

        if isfield(data, 'count')
            benchmark.count = data.count;
        else
            benchmark.count = [];
        end

        if isfield(data, 'limits')
            benchmark.limits = data.limits;
        else
            benchmark.limits = [];
        end

        if isfield(data, 'system')
            benchmark.systemInfo = data.system;
        else
            benchmark.systemInfo = [];
        end

        if isfield(data, 'fields')
            fields = data.fields;
        else
            error('JPL response for %s does not contain a "fields" entry.', c.tag);
        end

        if isstring(fields)
            fields = cellstr(fields);
        end

        if ischar(fields)
            fields = {fields};
        end

        fields = matlab.lang.makeValidName(fields);
        benchmark.fields = fields;

        if ~isfield(data, 'data')
            error('JPL response for %s does not contain a "data" entry.', c.tag);
        end

        rawRows = data.data;

        if isempty(rawRows)
            error('JPL response for %s returned zero rows.', c.tag);
        end

        % Convert JPL row data to numeric matrix
        %
        % The API usually returns rows as cell arrays of strings/numbers.
        % We convert everything to double row-by-row.
        nRows = numel(rawRows);
        nCols = numel(fields);
        A = nan(nRows, nCols);

        for i = 1:nRows
            row = rawRows{i};

            if iscell(row)
                rowNumeric = nan(1, numel(row));
                for j = 1:numel(row)
                    if isnumeric(row{j})
                        rowNumeric(j) = double(row{j});
                    else
                        rowNumeric(j) = str2double(string(row{j}));
                    end
                end
            elseif isnumeric(row)
                rowNumeric = double(row(:)).';
            else
                error('Unsupported row type in JPL response for %s.', c.tag);
            end

            if numel(rowNumeric) ~= nCols
                error('Row/field size mismatch in %s: row has %d cols, fields has %d.', ...
                    c.tag, numel(rowNumeric), nCols);
            end

            A(i, :) = rowNumeric;
        end

        benchmark.orbits = array2table(A, 'VariableNames', fields);
        benchmark.nOrbits = height(benchmark.orbits);

        % Some convenience summaries if those variables exist
        vars = benchmark.orbits.Properties.VariableNames;

        if any(strcmpi(vars, 'jacobi'))
            benchmark.jacobiMin = min(benchmark.orbits.jacobi);
            benchmark.jacobiMax = max(benchmark.orbits.jacobi);
        end

        if any(strcmpi(vars, 'period'))
            benchmark.periodMin = min(benchmark.orbits.period);
            benchmark.periodMax = max(benchmark.orbits.period);
        end

        if any(strcmpi(vars, 'stability'))
            benchmark.stabilityMin = min(benchmark.orbits.stability);
            benchmark.stabilityMax = max(benchmark.orbits.stability);
        end

        save(procMatFile, 'benchmark');

        fprintf('Saved raw JSON: %s\n', rawJsonFile);
        fprintf('Saved MAT file : %s\n', procMatFile);
        fprintf('Number of rows : %d\n', benchmark.nOrbits);
    end

    fprintf('\nAll selected JPL benchmark families were frozen successfully.\n');
end