clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/AstroToolbox')
startup

% ==========================================================
% Halo orbit differential correction demo
% Earth-Moon L1 northern halo family
% ==========================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('Halo Differential Correction Demo\n');
fprintf('============================================================\n');

% ----------------------------------------------------------
% Retrieve a JPL L1 northern halo seed
% ----------------------------------------------------------
data = astro.cr3bp.queryJPLPeriodicOrbits( ...
    'sys', 'earth-moon', ...
    'family', 'halo', ...
    'libr', 1, ...
    'branch', 'N');

seed = astro.cr3bp.parseJPLPeriodicOrbit(data, 1);
mu = seed.mu;

fprintf('\nJPL seed:\n');
fprintf('  x0     = %.12f\n', seed.state(1));
fprintf('  z0     = %.12f\n', seed.state(3));
fprintf('  vy0    = %.12f\n', seed.state(5));
fprintf('  Period = %.12f TU\n', seed.period);

% ----------------------------------------------------------
% Slightly perturb the seed before correction
% ----------------------------------------------------------
x0Guess = seed.state(1) + 1e-4;
z0Fix   = seed.state(3);          % use z0 as the family parameter
vy0Guess = seed.state(5) - 1e-4;
thalfGuess = seed.period / 2 * 1.0005;

% ----------------------------------------------------------
% Correct halo orbit
% ----------------------------------------------------------
out = astro.cr3bp.differentialCorrectionHalo( ...
    x0Guess, z0Fix, vy0Guess, thalfGuess, mu, 20, 1e-11);

fprintf('\nHalo correction result:\n');
fprintf('  Converged    : %d\n', out.converged);
fprintf('  Iterations   : %d\n', out.iterations);
fprintf('  Residual     : %.6e\n', out.residual);
fprintf('  x0           : %.12f\n', out.x0);
fprintf('  z0           : %.12f\n', out.z0);
fprintf('  vy0          : %.12f\n', out.vy0);
fprintf('  Half period  : %.12f TU\n', out.halfPeriod);
fprintf('  Full period  : %.12f TU\n', out.period);

C = astro.cr3bp.jacobiConstant(out.state0.', mu);
fprintf('  Corrected C  : %.12f\n', C);

% ----------------------------------------------------------
% Propagate corrected orbit over one full period
% ----------------------------------------------------------
opts.RelTol = 1e-12;
opts.AbsTol = 1e-12;
opts.Solver = 'ode113';

fullOut = astro.propagators.propagate( ...
    @(t,x) astro.cr3bp.eomCR3BP(t, x, mu), ...
    [0 out.period], out.state0, opts);

L = astro.cr3bp.lagrangePoints(mu);

% ----------------------------------------------------------
% 3D halo plot
% ----------------------------------------------------------
figure('Color','w');
hold on
grid on
box on
axis equal

plot3(fullOut.x(:,1), fullOut.x(:,2), fullOut.x(:,3), ...
    'b-', 'LineWidth', 1.8);

plot3(-mu, 0, 0, 'bo', 'MarkerSize', 10, 'LineWidth', 1.4, 'MarkerFaceColor', 'b');
plot3(1-mu, 0, 0, 'ko', 'MarkerSize', 7, 'LineWidth', 1.4, 'MarkerFaceColor', 'k');
plot3(L.L1(1), L.L1(2), 0, 'rs', 'MarkerSize', 7, 'LineWidth', 1.4);

xlabel('$x~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('$y~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
zlabel('$z~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('Corrected $L_1$ Northern Halo Orbit $(C = %.6f)$', C), ...
    'Interpreter', 'latex', 'FontSize', 15, 'FontWeight', 'bold');
legend('Halo orbit', 'Primary', 'Secondary', '$L_1$', ...
    'Interpreter', 'latex', 'Location', 'best');
view(35, 25)

% ----------------------------------------------------------
% 2D projections
% ----------------------------------------------------------
figure('Color','w');

subplot(1,2,1)
hold on
grid on
box on
axis equal
plot(fullOut.x(:,1), fullOut.x(:,2), 'b-', 'LineWidth', 1.5);
plot(-mu, 0, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
plot(1-mu, 0, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
plot(L.L1(1), L.L1(2), 'rs', 'MarkerSize', 6, 'LineWidth', 1.2);
xlabel('$x~[-]$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('$y~[-]$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
title('$x$-$y$ projection', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');

subplot(1,2,2)
hold on
grid on
box on
axis equal
plot(fullOut.x(:,1), fullOut.x(:,3), 'b-', 'LineWidth', 1.5);
plot(-mu, 0, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
plot(1-mu, 0, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
plot(L.L1(1), 0, 'rs', 'MarkerSize', 6, 'LineWidth', 1.2);
xlabel('$x~[-]$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('$z~[-]$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
title('$x$-$z$ projection', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');

fprintf('\nDone.\n');