function seed = getPlanarLyapunovSeed(sourceName, varargin)
%GETPLANARLYAPUNOVSEED Get a planar Lyapunov seed from JPL periodic orbit API.

    if nargin < 1 || isempty(sourceName)
        sourceName = 'jpl';
    end

    switch lower(string(sourceName))
        case "jpl"
            seed = localGetFromJPL(varargin{:});
        otherwise
            error('Unknown seed source: %s', sourceName);
    end
end

function seed = localGetFromJPL(varargin)
    p = inputParser;
    addParameter(p, 'sys', 'earth-moon', @(x)ischar(x)||isstring(x));
    addParameter(p, 'family', 'lyapunov', @(x)ischar(x)||isstring(x));
    addParameter(p, 'libr', 1, @(x)isnumeric(x) && isscalar(x));
    addParameter(p, 'branch', '', @(x)ischar(x)||isstring(x));
    addParameter(p, 'periodmin', [], @(x)isempty(x)||isnumeric(x));
    addParameter(p, 'periodmax', [], @(x)isempty(x)||isnumeric(x));
    addParameter(p, 'periodunits', '', @(x)ischar(x)||isstring(x));
    addParameter(p, 'jacobimin', [], @(x)isempty(x)||isnumeric(x));
    addParameter(p, 'jacobimax', [], @(x)isempty(x)||isnumeric(x));
    addParameter(p, 'stabmin', [], @(x)isempty(x)||isnumeric(x));
    addParameter(p, 'stabmax', [], @(x)isempty(x)||isnumeric(x));
    addParameter(p, 'index', 1, @(x)isnumeric(x) && isscalar(x) && x >= 1);
    parse(p, varargin{:});
    R = p.Results;

    data = astro.cr3bp.queryJPLPeriodicOrbits( ...
        'sys', R.sys, ...
        'family', R.family, ...
        'libr', R.libr, ...
        'branch', R.branch, ...
        'periodmin', R.periodmin, ...
        'periodmax', R.periodmax, ...
        'periodunits', R.periodunits, ...
        'jacobimin', R.jacobimin, ...
        'jacobimax', R.jacobimax, ...
        'stabmin', R.stabmin, ...
        'stabmax', R.stabmax);

    [fields, rows] = localExtractRows(data);

    idx = R.index;
    if idx > numel(rows)
        error('Requested row index %d, but only %d orbit rows were returned.', idx, numel(rows));
    end

    row = localNormalizeRow(rows{idx});

    seed = struct();
    seed.x0  = localGetFieldValue(row, fields, 'x');
    seed.y0  = localGetFieldValue(row, fields, 'y');
    seed.z0  = localGetFieldValue(row, fields, 'z');
    seed.vx0 = localGetFieldValue(row, fields, 'vx');
    seed.vy0 = localGetFieldValue(row, fields, 'vy');
    seed.vz0 = localGetFieldValue(row, fields, 'vz');
    seed.T   = localGetFieldValue(row, fields, 'period');
    seed.C   = localGetFieldValue(row, fields, 'jacobi');
    seed.stability = localGetFieldValue(row, fields, 'stability');
    seed.source = sprintf('JPL periodic orbit API (%s, %s, libr=%d, row=%d)', ...
        char(string(R.sys)), char(string(R.family)), R.libr, idx);
    seed.raw = data;

    if abs(seed.y0) > 1e-8 || abs(seed.z0) > 1e-8 || abs(seed.vx0) > 1e-8 || abs(seed.vz0) > 1e-8
        warning(['Returned JPL row is not exactly in planar symmetry form ', ...
                 '[x0;0;0;0;vy0;0].']);
    end
end

function [fields, rows] = localExtractRows(data)
    if ~isstruct(data)
        error('Unrecognized JPL response: expected a struct.');
    end

    if isfield(data,'fields') && isfield(data,'data')
        fields = string(data.fields);
        rows = data.data;
        return
    end

    error('Unrecognized JPL periodic orbit API response format.');
end

function row = localNormalizeRow(rawRow)
% Convert one JPL row into a flat cell array aligned with fields.

    row = rawRow;

    % Common case here: row is 1x1 cell, whose content is Nx1 cell
    while iscell(row) && numel(row) == 1
        row = row{1};
    end

    if ~iscell(row)
        error('JPL row could not be normalized to a cell array.');
    end

    % Force row orientation
    row = row(:).';
end

function val = localGetFieldValue(row, fields, targetName)
    idx = find(strcmpi(string(fields), targetName), 1);
    if isempty(idx)
        error('Field "%s" not found in JPL response.', targetName);
    end

    raw = row{idx};
    val = localToDouble(raw, targetName);
end

function val = localToDouble(raw, targetName)
    while iscell(raw) && numel(raw) == 1
        raw = raw{1};
    end

    if isnumeric(raw)
        if isempty(raw) || ~isscalar(raw)
            error('Field "%s" is numeric but not scalar.', targetName);
        end
        val = double(raw);
        return
    end

    if isstring(raw)
        if numel(raw) ~= 1
            error('Field "%s" is a string array, not scalar.', targetName);
        end
        val = str2double(raw);
        if isnan(val)
            error('Field "%s" string could not be converted to numeric.', targetName);
        end
        return
    end

    if ischar(raw)
        val = str2double(raw);
        if isnan(val)
            error('Field "%s" char value could not be converted to numeric.', targetName);
        end
        return
    end

    error('Unsupported field type for "%s": %s', targetName, class(raw));
end