function [seed1, seed2] = buildSeedFromSource(systemName, familyType, libr, branch, seedSource, opts)
%BUILDSEEDFROMSOURCE Build one or two corrected seeds for a periodic-orbit workflow.
%
% OUTPUTS
%   seed1, seed2 : structs with fields .state and .period
%
% Supported families:
%   - 'lyapunov'
%   - 'planar_lyapunov'
%   - 'halo'

if nargin < 6 || isempty(opts)
    opts = struct();
end

if ~isfield(opts,'validatedSeedName'); opts.validatedSeedName = ''; end
if ~isfield(opts,'jplIndex1'); opts.jplIndex1 = 1; end
if ~isfield(opts,'jplIndex2'); opts.jplIndex2 = 2; end

mu = astro.cr3bp.getSystemMu(systemName);

familyType = lower(string(familyType));
branch = lower(string(branch));

if familyType == "lyapunov"
    familyTypeNorm = "planar_lyapunov";
elseif familyType == "planar_lyapunov"
    familyTypeNorm = "planar_lyapunov";
elseif familyType == "halo"
    familyTypeNorm = "halo";
else
    error('Unsupported family type: %s', familyType);
end

switch familyTypeNorm
    case "planar_lyapunov"

        if strcmpi(seedSource, 'validated')
            if strlength(string(opts.validatedSeedName)) == 0
                error('validatedSeedName must be provided when seedSource = validated.');
            end

            S = astro.periodic.loadValidatedOrbit(opts.validatedSeedName, systemName);
            seed1.state = S.state0;
            seed1.period = S.T;

        elseif strcmpi(seedSource, 'jpl')
            data = astro.cr3bp.queryJPLPeriodicOrbits( ...
                'sys', strrep(systemName,'_','-'), ...
                'family', 'lyapunov', ...
                'libr', libr);

            fields = string(data.fields);
            rows = data.data;

            if opts.jplIndex1 > numel(rows) || opts.jplIndex2 > numel(rows)
                error('Requested JPL Lyapunov indices exceed available rows.');
            end

            seed1raw = localExtractSeed(rows{opts.jplIndex1}, fields);
            seed2raw = localExtractSeed(rows{opts.jplIndex2}, fields);

            corr1 = astro.cr3bp.differentialCorrectionPlanarLyapunov( ...
                seed1raw.x0, seed1raw.vy0, mu, 30, 1e-12);
            if ~corr1.converged
                error('Failed to correct first Lyapunov seed.');
            end

            corr2 = astro.cr3bp.differentialCorrectionPlanarLyapunov( ...
                seed2raw.x0, seed2raw.vy0, mu, 30, 1e-12);
            if ~corr2.converged
                error('Failed to correct second Lyapunov seed.');
            end

            seed1.state = corr1.state0;
            seed1.period = corr1.period;
            seed2.state = corr2.state0;
            seed2.period = corr2.period;
            return

        else
            error('Unsupported seedSource for Lyapunov family: %s', seedSource);
        end

        x0 = seed1.state(1);
        vy0 = seed1.state(5);
        thalf = 0.5 * seed1.period;

        uPred2 = [x0 + 1e-3; vy0; thalf];
        tangent0 = [1; 0; 0];
        tangent0 = tangent0 / norm(tangent0);

        corr2 = astro.cr3bp.correctPlanarLyapunovPseudoArc(uPred2, tangent0, mu, 20, 1e-11);
        if ~corr2.converged
            error('Failed to build second planar Lyapunov seed.');
        end

        seed2.state = corr2.state0;
        seed2.period = corr2.period;

    case "halo"

        if branch == "north"
            jplBranch = 'N';
        elseif branch == "south"
            jplBranch = 'S';
        else
            error('Halo branch must be north or south.');
        end

        data = astro.cr3bp.queryJPLPeriodicOrbits( ...
            'sys', strrep(systemName,'_','-'), ...
            'family', 'halo', ...
            'libr', libr, ...
            'branch', jplBranch);

        fields = string(data.fields);
        rows = data.data;

        if opts.jplIndex2 > numel(rows)
            error('Requested JPL halo indices exceed available rows.');
        end

        if strcmpi(seedSource, 'validated')
            if strlength(string(opts.validatedSeedName)) == 0
                error('validatedSeedName must be provided when seedSource = validated.');
            end

            S = astro.periodic.loadValidatedOrbit(opts.validatedSeedName, systemName);
            seed1.state = S.state0;
            seed1.period = S.T;

        elseif strcmpi(seedSource, 'jpl')
            if opts.jplIndex1 > numel(rows)
                error('Requested JPL halo indices exceed available rows.');
            end

            seed1raw = localExtractSeed(rows{opts.jplIndex1}, fields);
            corr1 = astro.cr3bp.differentialCorrectionHalo( ...
                seed1raw.x0, seed1raw.z0, seed1raw.vy0, 0.5*seed1raw.T, mu, 30, 1e-12);

            if ~corr1.converged
                error('Failed to correct first halo seed.');
            end

            seed1.state = corr1.state0;
            seed1.period = corr1.period;

        else
            error('Unsupported seedSource for halo family: %s', seedSource);
        end

        seed2raw = localExtractSeed(rows{opts.jplIndex2}, fields);
        corr2 = astro.cr3bp.differentialCorrectionHalo( ...
            seed2raw.x0, seed2raw.z0, seed2raw.vy0, 0.5*seed2raw.T, mu, 30, 1e-12);

        if ~corr2.converged
            error('Failed to correct second halo seed.');
        end

        seed2.state = corr2.state0;
        seed2.period = corr2.period;

    otherwise
        error('Unsupported family type: %s', familyType);
end
end

function seed = localExtractSeed(rawRow, fields)
    row = rawRow;
    while iscell(row) && numel(row) == 1
        row = row{1};
    end
    row = row(:).';

    seed = struct();
    seed.x0  = localGetFieldValue(row, fields, 'x');
    seed.y0  = localGetFieldValue(row, fields, 'y');
    seed.z0  = localGetFieldValue(row, fields, 'z');
    seed.vx0 = localGetFieldValue(row, fields, 'vx');
    seed.vy0 = localGetFieldValue(row, fields, 'vy');
    seed.vz0 = localGetFieldValue(row, fields, 'vz');
    seed.T   = localGetFieldValue(row, fields, 'period');
    seed.C   = localGetFieldValue(row, fields, 'jacobi');
end

function val = localGetFieldValue(row, fields, targetName)
    idx = find(strcmpi(fields, targetName), 1);
    raw = row{idx};

    while iscell(raw) && numel(raw) == 1
        raw = raw{1};
    end

    if isnumeric(raw)
        val = double(raw);
    elseif isstring(raw) || ischar(raw)
        val = str2double(raw);
    else
        error('Unsupported field type for "%s".', targetName);
    end
end