clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/AstroToolbox')
startup

% ==========================================================
% Stability analysis along a pseudo-arclength Lyapunov family
% Earth-Moon L1 planar Lyapunov family
% ==========================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('Family Stability Analysis Demo\n');
fprintf('============================================================\n');

% ----------------------------------------------------------
% Retrieve two nearby seeds from JPL
% ----------------------------------------------------------
data = astro.cr3bp.queryJPLPeriodicOrbits( ...
    'sys', 'earth-moon', ...
    'family', 'lyapunov', ...
    'libr', 1);

seed1 = astro.cr3bp.parseJPLPeriodicOrbit(data, 1);
seed2 = astro.cr3bp.parseJPLPeriodicOrbit(data, 2);

mu = seed1.mu;

% ----------------------------------------------------------
% Build family
% ----------------------------------------------------------
nMembers = 60;
ds = 5e-4;

family = astro.cr3bp.continueFamilyPseudoArc(seed1, seed2, nMembers, ds, mu);
nFam = numel(family);

if nFam < 2
    error('Too few family members generated for stability analysis.');
end

% ----------------------------------------------------------
% Analyze each family member
% ----------------------------------------------------------
Cvals = NaN(1,nFam);
Tvals = NaN(1,nFam);
x0vals = NaN(1,nFam);
vy0vals = NaN(1,nFam);

rhoVals = NaN(1,nFam);
logRhoVals = NaN(1,nFam);
nuVals = NaN(1,nFam);

lambdaUReal = NaN(1,nFam);
lambdaUImag = NaN(1,nFam);
lambdaSReal = NaN(1,nFam);
lambdaSImag = NaN(1,nFam);

center1Real = NaN(1,nFam);
center1Imag = NaN(1,nFam);
center2Real = NaN(1,nFam);
center2Imag = NaN(1,nFam);

for k = 1:nFam
    x0 = family(k).state0(:);
    T  = family(k).period;

    M = astro.cr3bp.monodromyMatrix(x0, T, mu);
    stab = astro.cr3bp.stabilityIndices(M);

    Cvals(k) = family(k).C;
    Tvals(k) = T;
    x0vals(k) = x0(1);
    vy0vals(k) = x0(5);

    rhoVals(k) = stab.spectralRadius;
    logRhoVals(k) = log10(stab.spectralRadius);
    nuVals(k) = real(stab.nu);

    lambdaUReal(k) = real(stab.unstableMultiplier);
    lambdaUImag(k) = imag(stab.unstableMultiplier);
    lambdaSReal(k) = real(stab.stableMultiplier);
    lambdaSImag(k) = imag(stab.stableMultiplier);

    if numel(stab.centerPair) >= 1
        center1Real(k) = real(stab.centerPair(1));
        center1Imag(k) = imag(stab.centerPair(1));
    end
    if numel(stab.centerPair) >= 2
        center2Real(k) = real(stab.centerPair(2));
        center2Imag(k) = imag(stab.centerPair(2));
    end
end

fprintf('\nFamily summary:\n');
fprintf('  Members            : %d\n', nFam);
fprintf('  C range            : [%.12f, %.12f]\n', min(Cvals), max(Cvals));
fprintf('  Period range [TU]  : [%.12f, %.12f]\n', min(Tvals), max(Tvals));
fprintf('  Spectral radius    : [%.12f, %.12f]\n', min(rhoVals), max(rhoVals));
fprintf('  Unstable mult.     : [%.12f, %.12f]\n', min(lambdaUReal), max(lambdaUReal));

% ----------------------------------------------------------
% Plot family in x-y plane
% ----------------------------------------------------------
L = astro.cr3bp.lagrangePoints(mu);

figure('Color','w');
hold on
axis equal
grid on
box on

cmap = parula(max(nFam,16));

for k = 1:nFam
    traj = family(k).traj;
    plot(traj(:,1), traj(:,2), 'Color', cmap(k,:), 'LineWidth', 1.2);
end

plot(-mu, 0, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
plot(1-mu, 0, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'k');
plot(L.L1(1), L.L1(2), 'rs', 'MarkerSize', 7, 'LineWidth', 1.4);
text(L.L1(1), L.L1(2), '  L1', 'FontSize', 10, 'FontWeight', 'bold');

colormap(cmap);
cb = colorbar;
ylabel(cb, 'Family index', 'FontSize', 12, 'FontWeight', 'bold');
clim([1 max(nFam,2)]);

xlabel('x [-]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('y [-]', 'FontSize', 13, 'FontWeight', 'bold');
title('L1 Lyapunov Family Used for Stability Analysis', ...
    'FontSize', 15, 'FontWeight', 'bold');

% ----------------------------------------------------------
% Period vs Jacobi
% ----------------------------------------------------------
figure('Color','w');
hold on
grid on
box on

plot(Cvals, Tvals, 'o-', 'LineWidth', 1.2, 'MarkerSize', 6);

xlabel('Jacobi constant C', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Period [TU]', 'FontSize', 13, 'FontWeight', 'bold');
title('Period vs Jacobi Constant Along Family', ...
    'FontSize', 15, 'FontWeight', 'bold');

% ----------------------------------------------------------
% Stability index nu vs family member
% ----------------------------------------------------------
figure('Color','w');
hold on
grid on
box on

plot(1:nFam, nuVals, 'o-', 'LineWidth', 1.2, 'MarkerSize', 6);

xlabel('Family member index', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Stability index \nu', 'FontSize', 13, 'FontWeight', 'bold');
title('Stability Index Along Family', ...
    'FontSize', 15, 'FontWeight', 'bold');

% ----------------------------------------------------------
% log10 spectral radius vs family member
% ----------------------------------------------------------
figure('Color','w');
hold on
grid on
box on

plot(1:nFam, logRhoVals, 'o-', 'LineWidth', 1.2, 'MarkerSize', 6);
yline(0.0, '--', 'LineWidth', 1.2);

xlabel('Family member index', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('log_{10}(\rho(M))', 'FontSize', 13, 'FontWeight', 'bold');
title('Log Spectral Radius Along Family', ...
    'FontSize', 15, 'FontWeight', 'bold');
legend('log_{10}(\rho)', 'Neutral threshold', 'Location', 'best');

% ----------------------------------------------------------
% Unstable / stable multipliers vs family member
% ----------------------------------------------------------
figure('Color','w');
hold on
grid on
box on

plot(1:nFam, lambdaUReal, 'o-', 'LineWidth', 1.2, 'MarkerSize', 6);
plot(1:nFam, lambdaSReal, 's-', 'LineWidth', 1.2, 'MarkerSize', 5);

xlabel('Family member index', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Real multiplier value', 'FontSize', 13, 'FontWeight', 'bold');
title('Unstable and Stable Reciprocal Multipliers Along Family', ...
    'FontSize', 15, 'FontWeight', 'bold');
legend('\lambda_u', '\lambda_s', 'Location', 'best');

% ----------------------------------------------------------
% Center pair in complex plane
% ----------------------------------------------------------
figure('Color','w');
hold on
grid on
box on
axis equal

th = linspace(0, 2*pi, 400);
plot(cos(th), sin(th), '--', 'LineWidth', 1.0);
plot(center1Real, center1Imag, 'o-', 'LineWidth', 1.2, 'MarkerSize', 6);
plot(center2Real, center2Imag, 's-', 'LineWidth', 1.0, 'MarkerSize', 5);

xlabel('Real(\lambda_c)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Imag(\lambda_c)', 'FontSize', 13, 'FontWeight', 'bold');
title('Center-Pair Eigenvalues Along Family', ...
    'FontSize', 15, 'FontWeight', 'bold');
legend('Unit circle', 'Center pair 1', 'Center pair 2', 'Location', 'best');

fprintf('\nDone.\n');