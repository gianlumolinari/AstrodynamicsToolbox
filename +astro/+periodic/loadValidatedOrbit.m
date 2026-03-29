function orbitData = loadValidatedOrbit(name, systemName)
%LOADVALIDATEDORBIT Load a validated periodic orbit from disk.
%
% USAGE
%   S = astro.periodic.loadValidatedOrbit(name)
%   S = astro.periodic.loadValidatedOrbit(name, systemName)

if nargin < 2
    systemName = '';
end

name = char(string(name));

if ~isempty(systemName)
    candidate = fullfile('data', 'validated', 'cr3bp', lower(systemName), [name '.mat']);
    if exist(candidate, 'file')
        orbitData = load(candidate);
        return
    end

    candidateValidated = fullfile('data', 'validated', 'cr3bp', lower(systemName), ...
        [name '_validated.mat']);
    if exist(candidateValidated, 'file')
        orbitData = load(candidateValidated);
        return
    end

    error('Validated orbit file not found for key "%s" in system "%s".', name, systemName);
end

baseDir = fullfile('data', 'validated', 'cr3bp');
if ~exist(baseDir, 'dir')
    error('Validated orbit database folder not found: %s', baseDir);
end

d = dir(baseDir);
d = d([d.isdir]);
d = d(~ismember({d.name}, {'.','..'}));

matches = {};

for k = 1:numel(d)
    folder = fullfile(baseDir, d(k).name);

    candidate = fullfile(folder, [name '.mat']);
    if exist(candidate, 'file')
        matches{end+1} = candidate; %#ok<AGROW>
    end

    candidateValidated = fullfile(folder, [name '_validated.mat']);
    if exist(candidateValidated, 'file')
        matches{end+1} = candidateValidated; %#ok<AGROW>
    end
end

if isempty(matches)
    error('Unknown validated orbit name: %s', name);
elseif numel(matches) > 1
    fprintf('Multiple matches found for "%s":\n', name);
    for k = 1:numel(matches)
        fprintf('  %d) %s\n', k, matches{k});
    end
    error('Orbit key "%s" is ambiguous. Call loadValidatedOrbit(name, systemName).', name);
else
    orbitData = load(matches{1});
end
end