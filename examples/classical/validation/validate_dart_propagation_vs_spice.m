clc;
clear;

cd('/Users/gianlucamolinari/Desktop/astroToolbox')
startup

% ==========================================================
% Minimal DART propagation vs SPICE validation
%
% Goal:
%   Compare toolbox Sun-centered two-body propagation against
%   DART SPICE truth over the available reconstructed interval.
% ==========================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('Minimal DART Propagation vs SPICE Validation\n');
fprintf('============================================================\n');

% ----------------------------------------------------------
% Paths
% ----------------------------------------------------------
repoRoot = '/Users/gianlucamolinari/Desktop/astroToolbox';
mkFile   = fullfile(repoRoot, 'data', 'spice', 'dart', 'spice_kernels', 'mk', 'dart_local.tm');

if exist(mkFile, 'file') ~= 2
    error('Meta-kernel not found: %s', mkFile);
end

% ----------------------------------------------------------
% Load kernels
% ----------------------------------------------------------
cspice_kclear
cspice_furnsh(mkFile)

fprintf('\nLoaded local DART meta-kernel successfully.\n');

dartID = cspice_bodn2c('DART');
fprintf('DART NAIF ID: %d\n', dartID);

% ----------------------------------------------------------
% SPK files available locally
% ----------------------------------------------------------
spkList = { ...
    fullfile(repoRoot, 'data', 'spice', 'dart', 'spice_kernels', 'spk', 'dart_2021_328_2022_091_rec_v01.bsp'), ...
    fullfile(repoRoot, 'data', 'spice', 'dart', 'spice_kernels', 'spk', 'dart_2022_091_2022_231_rec_v01.bsp'), ...
    fullfile(repoRoot, 'data', 'spice', 'dart', 'spice_kernels', 'spk', 'dart_2022_231_2022_269_rec_v03.bsp'), ...
    fullfile(repoRoot, 'data', 'spice', 'dart', 'spice_kernels', 'spk', 'dart_2022_269_2022_269_rec_v03.bsp'), ...
    fullfile(repoRoot, 'data', 'spice', 'dart', 'spice_kernels', 'spk', 'dart_2022_269_2022_269_spc_v04.bsp')};

spkList = spkList(cellfun(@(f) exist(f,'file') == 2, spkList));

if isempty(spkList)
    cspice_kclear
    error('No DART spacecraft SPK segments found.');
end

fprintf('Found %d DART spacecraft SPK segment(s).\n', numel(spkList));

% ----------------------------------------------------------
% Merge coverage windows
% ----------------------------------------------------------
room = 10000;
coverAll = [];

for iSeg = 1:numel(spkList)
    coverSeg = cspice_spkcov(spkList{iSeg}, dartID, room);
    if isempty(coverAll)
        coverAll = coverSeg;
    else
        coverAll = cspice_wnunid(coverAll, coverSeg);
    end
end

nIntervals = cspice_wncard(coverAll);
fprintf('Total merged coverage intervals: %d\n', nIntervals);

if nIntervals < 1
    cspice_kclear
    error('No valid coverage intervals found.');
end

[leftET, rightET] = cspice_wnfetd(coverAll, 1);

fprintf('Validation interval start: %s\n', cspice_et2utc(leftET,  'C', 3));
fprintf('Validation interval stop : %s\n', cspice_et2utc(rightET, 'C', 3));
fprintf('Validation interval span : %.6f days\n', (rightET-leftET)/86400);

% ----------------------------------------------------------
% Sample SPICE truth
% ----------------------------------------------------------
nSamples = 500;
etGrid = linspace(leftET, rightET, nSamples);

Xtruth = zeros(nSamples, 6);

for k = 1:nSamples
    stD = cspice_spkezr('DART', etGrid(k), 'J2000', 'NONE', 'SUN');
    Xtruth(k,:) = stD(:).';
end

x0 = Xtruth(1,:).';

% ----------------------------------------------------------
% Propagate with toolbox using Sun two-body dynamics
% ----------------------------------------------------------
muSun = astro.constants.SUN_MU;

odefun = @(t,x) twoBodySun(t, x, muSun);

opts.RelTol = 1e-12;
opts.AbsTol = 1e-12;
opts.Solver = 'ode113';

tSpanSec = etGrid - leftET;

out = astro.propagators.propagate(odefun, [tSpanSec(1), tSpanSec(end)], x0, opts);

% Interpolate to SPICE epochs
Xprop = zeros(nSamples, 6);
for j = 1:6
    Xprop(:,j) = interp1(out.t, out.x(:,j), tSpanSec, 'pchip');
end

% ----------------------------------------------------------
% Compute errors
% ----------------------------------------------------------
dr = Xprop(:,1:3) - Xtruth(:,1:3);
dv = Xprop(:,4:6) - Xtruth(:,4:6);

rErr = vecnorm(dr, 2, 2);
vErr = vecnorm(dv, 2, 2);

fprintf('\nError summary against SPICE truth:\n');
fprintf('  Initial position error [km]   : %.6e\n', rErr(1));
fprintf('  Final position error [km]     : %.6e\n', rErr(end));
fprintf('  Max position error [km]       : %.6e\n', max(rErr));
fprintf('  RMS position error [km]       : %.6e\n', sqrt(mean(rErr.^2)));

fprintf('  Initial velocity error [km/s] : %.6e\n', vErr(1));
fprintf('  Final velocity error [km/s]   : %.6e\n', vErr(end));
fprintf('  Max velocity error [km/s]     : %.6e\n', max(vErr));
fprintf('  RMS velocity error [km/s]     : %.6e\n', sqrt(mean(vErr.^2)));

fprintf('\nMinimal propagation-vs-SPICE validation completed successfully.\n');

cspice_kclear

% ==========================================================
% Local function: Sun-centered two-body dynamics
% ==========================================================
function dx = twoBodySun(~, x, muSun)
    r = x(1:3);
    v = x(4:6);

    rnorm = norm(r);
    a = -muSun * r / rnorm^3;

    dx = [v; a];
end