clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/AstroToolbox')
startup

% ==========================================================
% Halo family stability analysis demo
% Earth-Moon L1 northern halo family
% ==========================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('Halo Family Stability Analysis Demo\n');
fprintf('============================================================\n');

% ----------------------------------------------------------
% Retrieve two nearby JPL L1 northern halo seeds
% ----------------------------------------------------------
data = astro.cr3bp.queryJPLPeriodicOrbits( ...
    'sys', 'earth-moon', ...
    'family', 'halo', ...
    'libr', 1, ...
    'branch', 'N');

seed1 = astro.cr3bp.parseJPLPeriodicOrbit(data, 1);
seed2 = astro.cr3bp.parseJPLPeriodicOrbit(data, 2);

mu = seed1.mu;

% ----------------------------------------------------------
% Continue family
% ----------------------------------------------------------
nMembers = 30;
ds = 5e-4;

family = astro.cr3bp.continueHaloPseudoArc(seed1, seed2, nMembers, ds, mu);
nFam = numel(family);

if nFam < 2
    error('Too few halo family members generated.');
end

% ----------------------------------------------------------
% Stability analysis along family
% ----------------------------------------------------------
Cvals = NaN(1,nFam);
Tvals = NaN(1,nFam);
x0vals = NaN(1,nFam);
z0vals = NaN(1,nFam);
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
center1Abs  = NaN(1,nFam);
center2Abs  = NaN(1,nFam);
centerDrift = NaN(1,nFam);

for k = 1:nFam
    x0 = family(k).state0(:);
    T  = family(k).period;

    M = astro.cr3bp.monodromyMatrix(x0, T, mu);
    stab = astro.cr3bp.stabilityIndices(M);

    Cvals(k) = family(k).C;
    Tvals(k) = T;
    x0vals(k) = x0(1);
    z0vals(k) = x0(3);
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
        center1Abs(k)  = abs(stab.centerPair(1));
    end
    if numel(stab.centerPair) >= 2
        center2Real(k) = real(stab.centerPair(2));
        center2Imag(k) = imag(stab.centerPair(2));
        center2Abs(k)  = abs(stab.centerPair(2));
    end

    centerDrift(k) = max(abs([center1Abs(k)-1, center2Abs(k)-1]));
end

% ----------------------------------------------------------
% Finite-difference trend quantities
% ----------------------------------------------------------
dT_dC = gradient(Tvals, Cvals);
dz0_dC = gradient(z0vals, Cvals);
dlogrho_dC = gradient(logRhoVals, Cvals);
dlambdaU_dC = gradient(lambdaUReal, Cvals);

fprintf('\nHalo family summary:\n');
fprintf('  Stored members             : %d\n', nFam);
fprintf('  C range                    : [%.12f, %.12f]\n', min(Cvals), max(Cvals));
fprintf('  T range [TU]               : [%.12f, %.12f]\n', min(Tvals), max(Tvals));
fprintf('  z0 range                   : [%.12f, %.12f]\n', min(z0vals), max(z0vals));
fprintf('  Spectral radius range      : [%.12f, %.12f]\n', min(rhoVals), max(rhoVals));
fprintf('  log10(rho) range           : [%.12f, %.12f]\n', min(logRhoVals), max(logRhoVals));
fprintf('  Unstable multiplier range  : [%.12f, %.12f]\n', min(lambdaUReal), max(lambdaUReal));
fprintf('  Stable multiplier range    : [%.12f, %.12f]\n', min(lambdaSReal), max(lambdaSReal));
fprintf('  Max center-pair drift      : %.6e\n', max(centerDrift));

fprintf('\nDerived trend diagnostics:\n');
fprintf('  dT/dC range                : [%.6e, %.6e]\n', min(dT_dC), max(dT_dC));
fprintf('  dz0/dC range               : [%.6e, %.6e]\n', min(dz0_dC), max(dz0_dC));
fprintf('  dlog10(rho)/dC range       : [%.6e, %.6e]\n', min(dlogrho_dC), max(dlogrho_dC));
fprintf('  d(lambda_u)/dC range       : [%.6e, %.6e]\n', min(dlambdaU_dC), max(dlambdaU_dC));

L = astro.cr3bp.lagrangePoints(mu);
cmap = parula(max(nFam,16));

% ----------------------------------------------------------
% 3D family plot
% ----------------------------------------------------------
figure('Color','w');
hold on
grid on
box on
axis equal

for k = 1:nFam
    X = family(k).traj;
    plot3(X(:,1), X(:,2), X(:,3), 'Color', cmap(k,:), 'LineWidth', 1.2);
end

plot3(-mu, 0, 0, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
plot3(1-mu, 0, 0, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'k');
plot3(L.L1(1), 0, 0, 'rs', 'MarkerSize', 7, 'LineWidth', 1.4);

xlabel('$x~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('$y~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
zlabel('$z~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
title('$L_1$ Northern Halo Family Used for Stability Analysis', ...
    'Interpreter', 'latex', 'FontSize', 15, 'FontWeight', 'bold');
view(35, 25)

colormap(cmap);
cb = colorbar;
ylabel(cb, 'Family index', 'FontSize', 12, 'FontWeight', 'bold');
clim([1 max(nFam,2)]);

% ----------------------------------------------------------
% x-z projection
% ----------------------------------------------------------
figure('Color','w');
hold on
grid on
box on
axis equal

for k = 1:nFam
    X = family(k).traj;
    plot(X(:,1), X(:,3), 'Color', cmap(k,:), 'LineWidth', 1.2);
end

plot(-mu, 0, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
plot(1-mu, 0, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'k');
plot(L.L1(1), 0, 'rs', 'MarkerSize', 7, 'LineWidth', 1.4);

xlabel('$x~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('$z~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
title('$x$-$z$ Projection of the $L_1$ Halo Family', ...
    'Interpreter', 'latex', 'FontSize', 15, 'FontWeight', 'bold');

% ----------------------------------------------------------
% Stability dashboard
% ----------------------------------------------------------
figure('Color','w');

subplot(2,3,1)
hold on
grid on
box on
plot(Cvals, Tvals, 'o-', 'LineWidth', 1.2, 'MarkerSize', 5);
xlabel('$C$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('$T$ [TU]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
title('$T$ vs $C$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');

subplot(2,3,2)
hold on
grid on
box on
plot(Cvals, z0vals, 'o-', 'LineWidth', 1.2, 'MarkerSize', 5);
xlabel('$C$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('$z_0$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
title('$z_0$ vs $C$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');

subplot(2,3,3)
hold on
grid on
box on
plot(1:nFam, nuVals, 'o-', 'LineWidth', 1.2, 'MarkerSize', 5);
xlabel('Index', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('$\nu$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
title('Stability Index', 'Interpreter', 'latex','FontSize', 13, 'FontWeight', 'bold');

subplot(2,3,4)
hold on
grid on
box on
plot(1:nFam, logRhoVals, 'o-', 'LineWidth', 1.2, 'MarkerSize', 5);
yline(0.0, '--', 'LineWidth', 1.0);
xlabel('Index', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('$\log_{10}(\rho)$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
title('Log Spectral Radius', 'Interpreter', 'latex','FontSize', 13, 'FontWeight', 'bold');

subplot(2,3,5)
hold on
grid on
box on
plot(1:nFam, lambdaUReal, 'o-', 'LineWidth', 1.2, 'MarkerSize', 5);
plot(1:nFam, lambdaSReal, 's-', 'LineWidth', 1.1, 'MarkerSize', 4);
xlabel('Index', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Multiplier', 'FontSize', 12, 'FontWeight', 'bold');
title('$\lambda_u$ and $\lambda_s$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
legend('$\lambda_u$', '$\lambda_s$', 'Interpreter', 'latex', 'Location', 'best');

subplot(2,3,6)
hold on
grid on
box on
plot(1:nFam, centerDrift, 'o-', 'LineWidth', 1.2, 'MarkerSize', 5);
xlabel('Index', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('$\max(|\,|\lambda_c|-1\,|)$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
title('Center-Pair Drift from Unit Circle', 'Interpreter', 'latex','FontSize', 13, 'FontWeight', 'bold');

sgtitle('Halo Family Stability Dashboard','Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');

% ----------------------------------------------------------
% Derivative dashboard
% ----------------------------------------------------------
figure('Color','w');

subplot(2,2,1)
hold on
grid on
box on
plot(Cvals, dT_dC, 'o-', 'LineWidth', 1.2, 'MarkerSize', 5);
xlabel('$C$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('$dT/dC$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
title('Slope of Period Along Family', 'Interpreter', 'latex','FontSize', 13, 'FontWeight', 'bold');

subplot(2,2,2)
hold on
grid on
box on
plot(Cvals, dz0_dC, 'o-', 'LineWidth', 1.2, 'MarkerSize', 5);
xlabel('$C$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('$dz_0/dC$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
title('Slope of $z_0$ Along Family', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');

subplot(2,2,3)
hold on
grid on
box on
plot(Cvals, dlogrho_dC, 'o-', 'LineWidth', 1.2, 'MarkerSize', 5);
xlabel('$C$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('$d\log_{10}(\rho)/dC$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
title('Slope of Instability Measure', 'Interpreter', 'latex','FontSize', 13, 'FontWeight', 'bold');

subplot(2,2,4)
hold on
grid on
box on
plot(Cvals, dlambdaU_dC, 'o-', 'LineWidth', 1.2, 'MarkerSize', 5);
xlabel('$C$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('$d\lambda_u/dC$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
title('Slope of Unstable Multiplier', 'Interpreter', 'latex','FontSize', 13, 'FontWeight', 'bold');

sgtitle('Halo Family Derived Stability Trends', 'Interpreter', 'latex','FontSize', 16, 'FontWeight', 'bold');

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
plot(center1Real, center1Imag, 'o-', 'LineWidth', 1.2, 'MarkerSize', 5);
plot(center2Real, center2Imag, 's-', 'LineWidth', 1.0, 'MarkerSize', 4);

xlabel('Real$(\lambda_c)$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Imag$(\lambda_c)$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
title('Center-Pair Eigenvalues Along Halo Family', 'Interpreter', 'latex',...
    'FontSize', 15, 'FontWeight', 'bold');
legend('Unit circle', 'Center pair 1', 'Center pair 2', 'Location', 'best');

fprintf('\nDone.\n');