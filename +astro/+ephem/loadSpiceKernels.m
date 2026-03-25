function loadSpiceKernels(kernelDir)
%LOADSPICEKERNELS Loads SPICE kernels for planetary ephemerides or mission sets.
%
% INPUT
%   kernelDir : directory containing SPICE kernels
%
% Supported layouts:
%   1. Flat generic directory, e.g. data/spice
%   2. Nested mission directory, e.g. data/spice/juice
%
% Behaviour:
%   - clears previously loaded kernels
%   - searches recursively for relevant kernel files
%   - ensures at least one leap-seconds kernel is loaded
%   - loads planetary SPKs when present
%   - loads mission kernels when present
%
% Notes:
%   Call this once at the beginning of a script/session.

if nargin < 1
    kernelDir = fullfile(pwd, 'data', 'spice');
end

if ~exist(kernelDir, 'dir')
    error('SPICE kernel directory not found: %s', kernelDir);
end

cspice_kclear;

% -------------------------------------------------------------------------
% Recursively search all files
% -------------------------------------------------------------------------
allFiles = dir(fullfile(kernelDir, '**', '*'));
allFiles = allFiles(~[allFiles.isdir]);

allPaths = fullfile({allFiles.folder}, {allFiles.name});

% -------------------------------------------------------------------------
% Identify candidate kernels
% -------------------------------------------------------------------------
tlsIdx = endsWith(lower(allPaths), '.tls') | endsWith(lower(allPaths), '.tls.pc');
bspIdx = endsWith(lower(allPaths), '.bsp');
tfIdx  = endsWith(lower(allPaths), '.tf');
tpcIdx = endsWith(lower(allPaths), '.tpc');
bpcIdx = endsWith(lower(allPaths), '.bpc');
bcIdx  = endsWith(lower(allPaths), '.bc');
tscIdx = endsWith(lower(allPaths), '.tsc');

tlsFiles = allPaths(tlsIdx);
bspFiles = allPaths(bspIdx);
tfFiles  = allPaths(tfIdx);
tpcFiles = allPaths(tpcIdx);
bpcFiles = allPaths(bpcIdx);
bcFiles  = allPaths(bcIdx);
tscFiles = allPaths(tscIdx);

if isempty(tlsFiles)
    error('No leap-seconds kernel found in %s', kernelDir);
end

% Prefer generic naif0012.tls if present
tlsFile = '';
for k = 1:numel(tlsFiles)
    [~,name,ext] = fileparts(tlsFiles{k});
    if strcmpi([name ext], 'naif0012.tls')
        tlsFile = tlsFiles{k};
        break
    end
end
if isempty(tlsFile)
    tlsFile = tlsFiles{1};
end

% -------------------------------------------------------------------------
% Load kernels in a sensible order
% -------------------------------------------------------------------------
cspice_furnsh(tlsFile);

for k = 1:numel(tscFiles)
    cspice_furnsh(tscFiles{k});
end

for k = 1:numel(tpcFiles)
    cspice_furnsh(tpcFiles{k});
end

for k = 1:numel(bpcFiles)
    cspice_furnsh(bpcFiles{k});
end

for k = 1:numel(tfFiles)
    cspice_furnsh(tfFiles{k});
end

for k = 1:numel(bspFiles)
    cspice_furnsh(bspFiles{k});
end

for k = 1:numel(bcFiles)
    cspice_furnsh(bcFiles{k});
end

fprintf('Loaded SPICE kernels from:\n  %s\n', kernelDir);
fprintf('  TLS loaded : %s\n', tlsFile);
fprintf('  BSP count  : %d\n', numel(bspFiles));
fprintf('  FK count   : %d\n', numel(tfFiles));
fprintf('  PCK count  : %d\n', numel(tpcFiles) + numel(bpcFiles));
fprintf('  CK count   : %d\n', numel(bcFiles));
fprintf('  SCLK count : %d\n', numel(tscFiles));
end