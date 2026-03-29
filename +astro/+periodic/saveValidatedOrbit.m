function filePath = saveValidatedOrbit(state0, T, mu, meta)
%SAVEVALIDATEDORBIT Save a validated periodic orbit using a consistent naming scheme.
%
% INPUTS
%   state0 : 6x1 initial state
%   T      : period
%   mu     : CR3BP mass parameter
%   meta   : struct with fields
%       .system       e.g. 'earth_moon'
%       .family       e.g. 'halo', 'planar_lyapunov'
%       .libr         e.g. 1, 2
%       .branch       e.g. 'north', 'south' or ''
%       .qualifier    e.g. 'validated', 'small', 'medium', 'large'
%       .source       descriptive string
%       .closureError numeric
%
% OUTPUT
%   filePath : saved .mat file path

if nargin < 4 || isempty(meta)
    error('Metadata struct is required.');
end

requiredFields = {'system','family','libr','qualifier','source','closureError'};
for i = 1:numel(requiredFields)
    if ~isfield(meta, requiredFields{i})
        error('Missing meta.%s', requiredFields{i});
    end
end

if ~isfield(meta,'branch') || isempty(meta.branch)
    meta.branch = '';
end

baseDir = fullfile('data', 'validated', 'cr3bp', lower(meta.system));
if ~exist(baseDir, 'dir')
    mkdir(baseDir);
end

if isempty(meta.branch)
    fileName = sprintf('%s_L%d_%s.mat', ...
        lower(meta.family), meta.libr, lower(meta.qualifier));
else
    fileName = sprintf('%s_L%d_%s_%s.mat', ...
        lower(meta.family), meta.libr, lower(meta.branch), lower(meta.qualifier));
end

filePath = fullfile(baseDir, fileName);

family = meta.family;
libr = meta.libr;
branch = meta.branch;
source = meta.source;
closureError = meta.closureError;

save(filePath, 'state0', 'T', 'mu', 'family', 'libr', 'branch', 'source', 'closureError');

fprintf('Saved validated orbit to %s\n', filePath);
end