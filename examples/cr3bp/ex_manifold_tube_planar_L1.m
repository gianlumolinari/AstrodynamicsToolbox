clc;
clear;
close all;

% EX_MANIFOLD_TUBE_PLANAR_L1
% Test manifold generation from a known corrected planar Lyapunov orbit.

mu = 0.012150585609624;

% -------------------------------------------------------------------------
% KNOWN CORRECTED ORBIT
% Replace these with a trusted orbit from literature or from a previously
% validated computation.
% -------------------------------------------------------------------------
state0 = [4.0976123461511266E-1;
	-2.6988484146598425E-23;
	-2.9417515655701884E-26;
	-1.9237533891084223E-13;
	1.4666820372526499E+0;
	2.0898742783096624E-25];

T = 7.4458490878530990;   

if T <= 0
    error('Replace state0 and T with a known corrected planar Lyapunov orbit.');
end

orbit = astro.periodic.buildOrbitFromState(state0, T, mu);

fprintf('Using known corrected orbit:\n');
fprintf('  x0  = %.16f\n', orbit.state0(1));
fprintf('  vy0 = %.16f\n', orbit.state0(5));
fprintf('  T   = %.16f\n', orbit.T);

eps0 = 5e-5;
tfManifold = 12.0;

atlas = astro.manifolds.generateManifoldAtlas(orbit, mu, 'unstable', eps0, tfManifold, ...
    struct('numPhaseSamples',80, ...
           'includeBothSigns',true, ...
           'normalizeMode','position', ...
           'RelTol',1e-12, ...
           'AbsTol',1e-12));

atlas = astro.manifolds.sortTubeStrands(atlas, 'signphase');

fig = figure;
ax = axes(fig);
astro.plot.plotManifoldTube(atlas, struct('dimension','2d','ax',ax));
title('Planar L1 unstable manifold tube from known corrected orbit');

cerr = zeros(numel(atlas.branches),1);
for k = 1:numel(atlas.branches)
    cerr(k) = atlas.branches(k).CerrMax;
end

fprintf('Generated %d manifold branches.\n', numel(atlas.branches));
fprintf('Maximum Jacobi drift across all branches: %.3e\n', max(cerr));