function loadSpiceKernels(kernelDir)
%LOADSPICEKERNELS Loads required SPICE kernels for planetary ephemerides.
%
% INPUT
%   kernelDir : directory containing SPICE kernels
%
% Expected kernels:
%   naif0012.tls
%   de440s.bsp   (or de430.bsp / de440.bsp)
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

tlsFile = '';
bspFile = '';

tlsCandidates = {'naif0012.tls', 'naif0012.tls.pc'};
bspCandidates = {'de440.bsp', 'de430.bsp', 'de440s.bsp'};

for k = 1:numel(tlsCandidates)
    f = fullfile(kernelDir, tlsCandidates{k});
    if exist(f, 'file')
        tlsFile = f;
        break
    end
end

for k = 1:numel(bspCandidates)
    f = fullfile(kernelDir, bspCandidates{k});
    if exist(f, 'file')
        bspFile = f;
        break
    end
end

if isempty(tlsFile)
    error('No leap-seconds kernel found in %s', kernelDir);
end

if isempty(bspFile)
    error('No planetary SPK kernel found in %s', kernelDir);
end

cspice_furnsh(tlsFile);
cspice_furnsh(bspFile);

fprintf('Loaded SPICE kernels:\n');
fprintf('  TLS : %s\n', tlsFile);
fprintf('  BSP : %s\n', bspFile);
end