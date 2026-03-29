function familyData = loadValidatedFamily(systemName, familyName, libr, branch)
%LOADVALIDATEDFAMILY Load all validated family members matching a naming pattern.
%
% INPUTS
%   systemName : e.g. 'earth_moon'
%   familyName : e.g. 'halo', 'planar_lyapunov'
%   libr       : 1 or 2
%   branch     : optional, e.g. 'north', 'south'
%
% OUTPUT
%   familyData : struct array with loaded orbit files

if nargin < 4
    branch = '';
end

folder = fullfile('data', 'validated', 'cr3bp', lower(systemName));

if ~exist(folder, 'dir')
    error('Validated orbit folder does not exist: %s', folder);
end

if isempty(branch)
    pattern = sprintf('%s_L%d_*.mat', lower(familyName), libr);
else
    pattern = sprintf('%s_L%d_%s_*.mat', lower(familyName), libr, lower(branch));
end

files = dir(fullfile(folder, pattern));

familyData = struct([]);
for k = 1:numel(files)
    S = load(fullfile(folder, files(k).name));
    S.fileName = files(k).name;
    familyData(k,1) = S;
end

fprintf('\nLoaded %d validated family members matching %s\n', numel(familyData), pattern);
end