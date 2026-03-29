function files = listValidatedOrbits(systemName)
%LISTVALIDATEDORBITS List validated orbit files for a CR3BP system.
%
% INPUT
%   systemName : e.g. 'earth_moon'
%
% OUTPUT
%   files : struct array from dir()

if nargin < 1 || isempty(systemName)
    systemName = 'earth_moon';
end

folder = fullfile('data', 'validated', 'cr3bp', lower(systemName));

if ~exist(folder, 'dir')
    warning('Validated orbit folder does not exist: %s', folder);
    files = struct([]);
    return
end

files = dir(fullfile(folder, '*.mat'));

fprintf('\nValidated orbit files in %s\n', folder);
fprintf('----------------------------------------\n');
for k = 1:numel(files)
    fprintf('%2d) %s\n', k, files(k).name);
end
end