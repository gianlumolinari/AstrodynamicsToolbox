clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/AstroToolbox')
startup

fprintf('\n============================================================\n');
fprintf('Howell 1984 Halo Validation (Table I, L1 family, mu = 0.04)\n');
fprintf('============================================================\n');

% ==========================================================
% Howell 1984, Table I: Initial conditions for L1 family at mu = 0.04
%
% State convention:
%   X0 = [x0; 0; z0; 0; yd0; 0]
% ==========================================================

mu = 0.04;

x0_tab = [ ...
    0.723268, ...
    0.729988, ...
    0.753700, ...
    0.777413, ...
    0.801125, ...
    0.817724];

z0_tab = [ ...
    0.040000, ...
    0.215589, ...
    0.267595, ...
    0.284268, ...
    0.299382, ...
    0.313788];

yd0_tab = [ ...
    0.198019, ...
    0.397259, ...
    0.399909, ...
    0.361870, ...
    0.312474, ...
    0.271306];

thalf_tab = [ ...
    1.300177, ...
    1.348532, ...
    1.211253, ...
    1.101099, ...
    1.017241, ...
    0.978635];

C_tab = [ ...
    3.329168, ...
    3.030033, ...
    2.937178, ...
    2.928754, ...
    2.930700, ...
    2.929481];

nCases = numel(x0_tab);

% Storage
C_comp = NaN(1,nCases);
C_err = NaN(1,nCases);
period_tab = 2*thalf_tab;
closure_pos = NaN(1,nCases);
closure_vel = NaN(1,nCases);

corr_x0_err   = NaN(1,nCases);
corr_z0_err   = NaN(1,nCases);
corr_yd0_err  = NaN(1,nCases);
corr_th_err   = NaN(1,nCases);
corr_C_err    = NaN(1,nCases);
corr_residual = NaN(1,nCases);
corr_conv     = false(1,nCases);

opts.RelTol = 1e-12;
opts.AbsTol = 1e-12;
opts.Solver = 'ode113';

trajCell = cell(1,nCases);

fprintf('\nTabulated-case checks:\n');

for k = 1:nCases
    Xref = [x0_tab(k); 0; z0_tab(k); 0; yd0_tab(k); 0];
    Tref = 2*thalf_tab(k);

    % ------------------------------------------------------
    % 1. Jacobi constant from tabulated IC
    % ------------------------------------------------------
    C_comp(k) = astro.cr3bp.jacobiConstant(Xref.', mu);
    C_err(k) = C_comp(k) - C_tab(k);

    % ------------------------------------------------------
    % 2. Propagate one full period and measure closure
    % ------------------------------------------------------
    out = astro.propagators.propagate( ...
        @(t,x) astro.cr3bp.eomCR3BP(t, x, mu), ...
        [0 Tref], Xref, opts);

    Xend = out.x(end,:).';
    closure_pos(k) = norm(Xend(1:3) - Xref(1:3));
    closure_vel(k) = norm(Xend(4:6) - Xref(4:6));

    trajCell{k} = out.x;

    fprintf('\nCase %d\n', k);
    fprintf('  Table C              = %.12f\n', C_tab(k));
    fprintf('  Computed C           = %.12f\n', C_comp(k));
    fprintf('  C error              = %.3e\n', C_err(k));
    fprintf('  Full period          = %.12f TU\n', Tref);
    fprintf('  Closure pos error    = %.3e\n', closure_pos(k));
    fprintf('  Closure vel error    = %.3e\n', closure_vel(k));

    % ------------------------------------------------------
    % 3. Differential-correction recovery from perturbed guess
    % ------------------------------------------------------
    x0_guess    = x0_tab(k) * (1 + 1e-4);
    z0_fixed    = z0_tab(k);
    yd0_guess   = yd0_tab(k) * (1 - 1e-4);
    thalf_guess = thalf_tab(k) * (1 + 5e-4);

    corr = astro.cr3bp.differentialCorrectionHalo( ...
        x0_guess, z0_fixed, yd0_guess, thalf_guess, mu, 25, 1e-11);

    corr_conv(k) = corr.converged;
    corr_residual(k) = corr.residual;

    Xcorr = corr.state0(:);
    Ccorr = astro.cr3bp.jacobiConstant(Xcorr.', mu);

    corr_x0_err(k)  = Xcorr(1) - x0_tab(k);
    corr_z0_err(k)  = Xcorr(3) - z0_tab(k);
    corr_yd0_err(k) = Xcorr(5) - yd0_tab(k);
    corr_th_err(k)  = corr.halfPeriod - thalf_tab(k);
    corr_C_err(k)   = Ccorr - C_tab(k);

    fprintf('  Corrector converged  = %d\n', corr.converged);
    fprintf('  Corrector residual   = %.3e\n', corr.residual);
    fprintf('  x0 error             = %.3e\n', corr_x0_err(k));
    fprintf('  z0 error             = %.3e\n', corr_z0_err(k));
    fprintf('  yd0 error            = %.3e\n', corr_yd0_err(k));
    fprintf('  T/2 error            = %.3e\n', corr_th_err(k));
    fprintf('  Corrected C error    = %.3e\n', corr_C_err(k));
end

% ==========================================================
% Summary table to command window
% ==========================================================
fprintf('\n============================================================\n');
fprintf('Summary\n');
fprintf('============================================================\n');
fprintf('Max |C_table - C_computed|      = %.3e\n', max(abs(C_err)));
fprintf('Max closure position error      = %.3e\n', max(closure_pos));
fprintf('Max closure velocity error      = %.3e\n', max(closure_vel));
fprintf('Max |x0 correction error|       = %.3e\n', max(abs(corr_x0_err)));
fprintf('Max |z0 correction error|       = %.3e\n', max(abs(corr_z0_err)));
fprintf('Max |yd0 correction error|      = %.3e\n', max(abs(corr_yd0_err)));
fprintf('Max |T/2 correction error|      = %.3e\n', max(abs(corr_th_err)));
fprintf('Max |corrected C error|         = %.3e\n', max(abs(corr_C_err)));
fprintf('All correctors converged        = %d\n', all(corr_conv));

% ==========================================================
% Plot: 3D reference family
% ==========================================================
L = astro.cr3bp.lagrangePoints(mu);
cmap = parula(max(nCases,16));

figure('Color','w');
hold on
grid on
box on
axis equal

for k = 1:nCases
    X = trajCell{k};
    plot3(X(:,1), X(:,2), X(:,3), 'Color', cmap(k,:), 'LineWidth', 1.4);
end

plot3(-mu, 0, 0, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
plot3(1-mu, 0, 0, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'k');
plot3(L.L1(1), 0, 0, 'rs', 'MarkerSize', 7, 'LineWidth', 1.4);

xlabel('$x~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('$y~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
zlabel('$z~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
title('Howell 1984 Table I: $L_1$ Halo Family at $\mu=0.04$', ...
    'Interpreter', 'latex', 'FontSize', 15, 'FontWeight', 'bold');
view(35,25)

colormap(cmap);
cb = colorbar;
ylabel(cb, 'Table column index', 'FontSize', 12, 'FontWeight', 'bold');
clim([1 nCases]);

% ==========================================================
% Plot: x-z projection
% ==========================================================
figure('Color','w');
hold on
grid on
box on
axis equal

for k = 1:nCases
    X = trajCell{k};
    plot(X(:,1), X(:,3), 'Color', cmap(k,:), 'LineWidth', 1.4);
end

plot(-mu, 0, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
plot(1-mu, 0, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'k');
plot(L.L1(1), 0, 'rs', 'MarkerSize', 7, 'LineWidth', 1.4);

xlabel('$x~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('$z~[-]$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
title('Howell 1984 Table I: $x$-$z$ Projection', ...
    'Interpreter', 'latex', 'FontSize', 15, 'FontWeight', 'bold');

% ==========================================================
% Plot: validation diagnostics dashboard
% ==========================================================
figure('Color','w');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

nexttile
hold on
grid on
box on
plot(1:nCases, abs(C_err), 'o-', 'LineWidth', 1.2, 'MarkerSize', 6);
xlabel('Table column', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('$|C_{\mathrm{comp}}-C_{\mathrm{tab}}|$', ...
    'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
title('Jacobi Constant Consistency', 'FontSize', 13, 'FontWeight', 'bold');

nexttile
hold on
grid on
box on
plot(1:nCases, closure_pos, 'o-', 'LineWidth', 1.2, 'MarkerSize', 6);
plot(1:nCases, closure_vel, 's-', 'LineWidth', 1.2, 'MarkerSize', 6);
xlabel('Table column', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Error norm', 'FontSize', 12, 'FontWeight', 'bold');
title('One-Period Orbit Closure', 'FontSize', 13, 'FontWeight', 'bold');
legend('Position closure', 'Velocity closure', 'Location', 'best');

nexttile
hold on
grid on
box on
plot(1:nCases, abs(corr_x0_err), 'o-', 'LineWidth', 1.2, 'MarkerSize', 6);
plot(1:nCases, abs(corr_z0_err), 's-', 'LineWidth', 1.2, 'MarkerSize', 6);
plot(1:nCases, abs(corr_yd0_err), 'd-', 'LineWidth', 1.2, 'MarkerSize', 6);
xlabel('Table column', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Absolute error', 'FontSize', 12, 'FontWeight', 'bold');
title('Recovery of Tabulated Initial Conditions', 'FontSize', 13, 'FontWeight', 'bold');
legend('$x_0$', '$z_0$', '$\dot y_0$', 'Interpreter', 'latex', 'Location', 'best');

nexttile
hold on
grid on
box on
plot(1:nCases, abs(corr_th_err), 'o-', 'LineWidth', 1.2, 'MarkerSize', 6);
plot(1:nCases, abs(corr_C_err), 's-', 'LineWidth', 1.2, 'MarkerSize', 6);
plot(1:nCases, corr_residual, 'd-', 'LineWidth', 1.2, 'MarkerSize', 6);
xlabel('Table column', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Magnitude', 'FontSize', 12, 'FontWeight', 'bold');
title('Corrector Accuracy', 'FontSize', 13, 'FontWeight', 'bold');
legend('$|T/2-T/2_{\mathrm{tab}}|$', '$|C-C_{\mathrm{tab}}|$', 'Residual', ...
    'Interpreter', 'latex', 'Location', 'best');

fprintf('\nDone.\n');