clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/AstroToolbox')
startup

% ==========================================================
% Halo eigenspectrum inspection demo
% Earth-Moon L1 northern halo family
% ==========================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('Halo Eigenspectrum Inspection Demo\n');
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

if nFam < 3
    error('Too few halo family members generated.');
end

% Select representative members
idxList = unique([1, round(nFam/2), nFam]);
nSel = numel(idxList);

fprintf('\nInspecting halo family members:\n');
disp(idxList);

eigStore = cell(1,nSel);
stabStore = cell(1,nSel);

for j = 1:nSel
    k = idxList(j);

    x0 = family(k).state0(:);
    T  = family(k).period;

    M = astro.cr3bp.monodromyMatrix(x0, T, mu);
    stab = astro.cr3bp.stabilityIndices(M);

    eigStore{j} = stab.eigvals;
    stabStore{j} = stab;

    fprintf('\nMember %d\n', k);
    fprintf('  C                  = %.12f\n', family(k).C);
    fprintf('  T                  = %.12f TU\n', family(k).period);
    fprintf('  Spectral radius    = %.12f\n', stab.spectralRadius);
    fprintf('  Unstable multiplier= %.12f%+.12fi\n', ...
        real(stab.unstableMultiplier), imag(stab.unstableMultiplier));
    fprintf('  Stable multiplier  = %.12f%+.12fi\n', ...
        real(stab.stableMultiplier), imag(stab.stableMultiplier));
    fprintf('  nu                 = %.12f%+.12fi\n', real(stab.nu), imag(stab.nu));

    fprintf('  Full eigenspectrum:\n');
    for m = 1:numel(stab.eigvals)
        lam = stab.eigvals(m);
        fprintf('    lambda_%d = %.12f%+.12fi   |lambda| = %.12f\n', ...
            m, real(lam), imag(lam), abs(lam));
    end
end

% ----------------------------------------------------------
% Plot eigenspectra in complex plane
% ----------------------------------------------------------
figure('Color','w');

th = linspace(0, 2*pi, 400);
ucx = cos(th);
ucy = sin(th);

for j = 1:nSel
    subplot(1,nSel,j)
    hold on
    grid on
    box on
    axis equal

    plot(ucx, ucy, 'k--', 'LineWidth', 1.0);

    lam = eigStore{j};
    stab = stabStore{j};

    plot(real(lam), imag(lam), 'ko', 'MarkerSize', 7, 'LineWidth', 1.2);

    plot(real(stab.unstableMultiplier), imag(stab.unstableMultiplier), ...
        'ro', 'MarkerSize', 8, 'LineWidth', 1.6);
    plot(real(stab.stableMultiplier), imag(stab.stableMultiplier), ...
        'bs', 'MarkerSize', 8, 'LineWidth', 1.6);

    if all(~isnan(stab.centerPair))
        plot(real(stab.centerPair), imag(stab.centerPair), ...
            'gd', 'MarkerSize', 8, 'LineWidth', 1.6);
    end

    title(sprintf('Member %d', idxList(j)), 'FontSize', 13, 'FontWeight', 'bold');
    xlabel('Re$(\lambda)$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Im$(\lambda)$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
end

sgtitle('Halo Monodromy Eigenspectra', 'FontSize', 16, 'FontWeight', 'bold');

% ----------------------------------------------------------
% Plot eigenvalue moduli by member
% ----------------------------------------------------------
eigAbsMat = NaN(6, nSel);
for j = 1:nSel
    eigAbsMat(:,j) = sort(abs(eigStore{j}), 'descend');
end

figure('Color','w');
hold on
grid on
box on

for i = 1:6
    plot(idxList, eigAbsMat(i,:), 'o-', 'LineWidth', 1.2, 'MarkerSize', 6);
end

yline(1.0, 'k--', 'LineWidth', 1.0);

xlabel('Family member index', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('$|\lambda|$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
title('Sorted Eigenvalue Moduli for Selected Halo Members', ...
    'FontSize', 15, 'FontWeight', 'bold');

fprintf('\nDone.\n');