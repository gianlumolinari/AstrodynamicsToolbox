clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/astroToolbox')
startup

mkFile  = '/Users/gianlucamolinari/Desktop/astroToolbox/data/spice/dart/spice_kernels/mk/dart_local.tm';
dartSpk = '/Users/gianlucamolinari/Desktop/astroToolbox/data/spice/dart/spice_kernels/spk/dart_2022_269_2022_269_spc_v04.bsp';

cspice_kclear
cspice_furnsh(mkFile)

fprintf('\nLoaded DART local meta-kernel successfully.\n');

dartID = cspice_bodn2c('DART');
fprintf('DART NAIF ID: %d\n', dartID);

room = 1000;
cover = cspice_spkcov(dartSpk, dartID, room);
nIntervals = cspice_wncard(cover);

fprintf('Coverage intervals found: %d\n', nIntervals);

if nIntervals < 1
    cspice_kclear
    error('No DART coverage found in spacecraft SPK.');
end

[leftET, rightET] = cspice_wnfetd(cover, 1);

fprintf('Coverage start : %s\n', cspice_et2utc(leftET,  'C', 3));
fprintf('Coverage stop  : %s\n', cspice_et2utc(rightET, 'C', 3));
fprintf('Coverage hours : %.3f\n', (rightET-leftET)/3600);

nSamples = 300;
etGrid = linspace(leftET, rightET, nSamples);

rDart  = zeros(nSamples,3);
rEarth = zeros(nSamples,3);

for k = 1:nSamples
    et = etGrid(k);

    stD = cspice_spkezr('DART', et, 'J2000', 'NONE', 'SUN');
    stE = cspice_spkezr('EARTH', et, 'J2000', 'NONE', 'SUN');

    rDart(k,:)  = stD(1:3).';
    rEarth(k,:) = stE(1:3).';
end



cspice_kclear