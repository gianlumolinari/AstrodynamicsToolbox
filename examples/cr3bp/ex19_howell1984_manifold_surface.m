clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/AstroToolbox')
startup

fprintf('\n============================================================\n');
fprintf('Howell 1984 Manifold Surface Example\n');
fprintf('============================================================\n');

% ==========================================================
% Howell 1984 Table I, L1 family at mu = 0.04
% Choose one validated benchmark case
% ==========================================================
mu = 0.04;

% Recommended starting case: column 2
x0_tab    = 0.729988;
z0_tab    = 0.215589;
yd0_tab   = 0.397259;
thalf_tab = 1.348532;
C_tab     = 3.030033;

x0 = [x0_tab; 0; z0_tab; 0; yd0_tab; 0];
T  = 2*thalf_tab;

fprintf('\nSelected Howell benchmark case:\n');
fprintf('  x0     = %.12f\n', x0_tab);
fprintf('  z0     = %.12f\n', z0_tab);
fprintf('  yd0    = %.12f\n', yd0_tab);
fprintf('  T/2    = %.12f\n', thalf_tab);
fprintf('  T      = %.12f\n', T);
fprintf('  C(tab) = %.12f\n', C_tab);

Ccomp = astro.cr3bp.jacobiConstant(x0.', mu);
fprintf('  C(comp)= %.12f\n', Ccomp);
fprintf('  |dC|   = %.3e\n', abs(Ccomp - C_tab));

% ----------------------------------------------------------
% Reference periodic orbit
% ----------------------------------------------------------
opts.RelTol = 1e-12;
opts.AbsTol = 1e-12;
opts.Solver = 'ode113';

orbit = astro.propagators.propagate( ...
    @(t,x) astro.cr3bp.eomCR3BP(t, x, mu), ...
    [0 T], x0, opts);

% ----------------------------------------------------------
% Manifold settings
% ----------------------------------------------------------
branchType = 'stable';   % 'stable' or 'unstable'
sense = -1;              % +1 or -1
nPhase = 120;
epsManifold = 1e-6;

Xseed = astro.cr3bp.howellManifoldBranch( ...
    x0, T, mu, nPhase, branchType, sense, epsManifold);

switch lower(branchType)
    case 'stable'
        tMan = -2.0*T;
        manColor = [0.10 0.25 0.95];
    case 'unstable'
        tMan =  2.0*T;
        manColor = [0.85 0.10 0.10];
end

traj = cell(nPhase,1);
fprintf('\nPropagating manifold surface...\n');
for k = 1:nPhase
    traj{k} = astro.cr3bp.propagateManifold(Xseed(k,:).', tMan, mu);
end

L = astro.cr3bp.lagrangePoints(mu);

% ----------------------------------------------------------
% Figure 1: Ross/Howell-style x-y projection
% ----------------------------------------------------------
figure('Color','w');
hold on
grid on
box on
axis equal

for k = 1:nPhase
    X = traj{k}.x;
    plot(X(:,1), X(:,2), '-', 'Color', manColor, 'LineWidth', 0.55);
end

plot(orbit.x(:,1), orbit.x(:,2), 'r-', 'LineWidth', 1.8);

plot(-mu, 0, 'o', 'MarkerSize', 7, ...
    'MarkerFaceColor', [0.55 0.35 0.15], ...
    'MarkerEdgeColor', [0.55 0.35 0.15]);

plot(1-mu, 0, 'o', 'MarkerSize', 8, ...
    'MarkerFaceColor', [0.75 0.75 0.75], ...
    'MarkerEdgeColor', [0.75 0.75 0.75]);

xlabel('X (DU)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Y (DU)', 'FontSize', 16, 'FontWeight', 'bold');
title(sprintf('Howell 1984 %s Manifold : Sense %d', ...
    capitalizeFirst(branchType), sense), ...
    'FontSize', 18, 'FontWeight', 'bold');

text(-0.05, 0.08, 'EARTH', 'FontSize', 16, 'FontWeight', 'bold');
text(0.95, -0.10, 'MOON', 'FontSize', 16, 'FontWeight', 'bold');
text(0.42, 0.42, 'HALO ORBIT', 'FontSize', 16, 'FontWeight', 'bold');

% ----------------------------------------------------------
% Figure 2: x-z projection
% ----------------------------------------------------------
figure('Color','w');
hold on
grid on
box on
axis equal

for k = 1:nPhase
    X = traj{k}.x;
    plot(X(:,1), X(:,3), '-', 'Color', manColor, 'LineWidth', 0.55);
end

plot(orbit.x(:,1), orbit.x(:,3), 'r-', 'LineWidth', 1.8);

plot(-mu, 0, 'o', 'MarkerSize', 7, ...
    'MarkerFaceColor', [0.55 0.35 0.15], ...
    'MarkerEdgeColor', [0.55 0.35 0.15]);

plot(1-mu, 0, 'o', 'MarkerSize', 8, ...
    'MarkerFaceColor', [0.75 0.75 0.75], ...
    'MarkerEdgeColor', [0.75 0.75 0.75]);

xlabel('X (DU)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Z (DU)', 'FontSize', 16, 'FontWeight', 'bold');
title(sprintf('Howell 1984 %s Manifold in X-Z Projection', ...
    capitalizeFirst(branchType)), ...
    'FontSize', 18, 'FontWeight', 'bold');

fprintf('\nDone.\n');

function s = capitalizeFirst(str)
str = char(string(str));
if isempty(str)
    s = str;
else
    s = [upper(str(1)) str(2:end)];
end
end