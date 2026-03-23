clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/AstroToolbox')
startup

% ==========================================================
% JPL full-family retrieval and propagation demo
%
% Retrieves a CR3BP family from the JPL periodic-orbits API,
% propagates every returned orbit with our own CR3BP equations,
% and plots the full family using orbit index as color variable.
%
% Default example:
%   Earth-Moon L1 northern halo family
% ==========================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('JPL Full Family Retrieval and Propagation Demo\n');
fprintf('============================================================\n');

% ----------------------------------------------------------
% User settings
% ----------------------------------------------------------
sysName = 'earth-moon';
familyName = 'lyapunov';
librPoint = 1;
branchName = '';

% Optional filter
periodMaxTU = 3.0;

% Downsample family if very large
maxOrbitsToPlot = inf;   % e.g. 50 if needed

% Propagation options
opts.RelTol = 1e-12;
opts.AbsTol = 1e-12;
opts.Solver = 'ode113';

% ----------------------------------------------------------
% Query JPL family
% ----------------------------------------------------------
data = astro.cr3bp.queryJPLPeriodicOrbits( ...
    'sys', sysName, ...
    'family', familyName, ...
    'libr', librPoint, ...
    'branch', branchName, ...
    'periodmax', periodMaxTU, ...
    'periodunits', 'TU');

nTotal = str2double(data.count);
fprintf('\nJPL API returned %d family members.\n', nTotal);

nUse = min(nTotal, maxOrbitsToPlot);
fprintf('Using %d orbits for propagation/plotting.\n', nUse);

% ----------------------------------------------------------
% Parse all selected family members
% ----------------------------------------------------------
orbits = repmat(struct( ...
    'state', zeros(6,1), ...
    'jacobi', NaN, ...
    'period', NaN, ...
    'stability', NaN, ...
    'mu', NaN, ...
    'lunit_km', NaN, ...
    'tunit_s', NaN, ...
    'systemName', '', ...
    'family', '', ...
    'branch', '', ...
    'librationPoint', ''), nUse, 1);

for k = 1:nUse
    orbits(k) = astro.cr3bp.parseJPLPeriodicOrbit(data, k);
end

mu = orbits(1).mu;
L = astro.cr3bp.lagrangePoints(mu);

% ----------------------------------------------------------
% Orbit-index coloring
% ----------------------------------------------------------
cvals = 1:nUse;
cLabel = 'Orbit index';

cmap = parula(max(nUse, 16));
cmin = 1;
cmax = max(nUse, 2);

fprintf('\nFamily summary:\n');
fprintf('  System           : %s\n', orbits(1).systemName);
fprintf('  Family           : %s\n', orbits(1).family);
fprintf('  Branch           : %s\n', orbits(1).branch);
fprintf('  Libration point  : %s\n', orbits(1).librationPoint);
fprintf('  mu               : %.15f\n', mu);
fprintf('  Number of orbits : %d\n', nUse);

% ----------------------------------------------------------
% Plot full family in x-y plane
% ----------------------------------------------------------
figure('Color','w');
hold on
axis equal
grid on
box on

plot(-mu, 0, 'bo', 'MarkerSize', 10, 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
plot(1-mu, 0, 'ko', 'MarkerSize', 7, 'LineWidth', 1.5, 'MarkerFaceColor', 'k');

plot(L.L1(1), L.L1(2), 'rs', 'MarkerSize', 7, 'LineWidth', 1.4);
plot(L.L2(1), L.L2(2), 'rs', 'MarkerSize', 7, 'LineWidth', 1.4);
plot(L.L3(1), L.L3(2), 'rs', 'MarkerSize', 7, 'LineWidth', 1.4);
plot(L.L4(1), L.L4(2), 'ms', 'MarkerSize', 7, 'LineWidth', 1.4);
plot(L.L5(1), L.L5(2), 'ms', 'MarkerSize', 7, 'LineWidth', 1.4);

text(L.L1(1), L.L1(2), '  L1', 'FontSize', 10, 'FontWeight', 'bold');
text(L.L2(1), L.L2(2), '  L2', 'FontSize', 10, 'FontWeight', 'bold');

for k = 1:nUse
    x0 = orbits(k).state;
    T = orbits(k).period;

    out = astro.propagators.propagate( ...
        @(t,x) astro.cr3bp.eomCR3BP(t, x, mu), ...
        [0, T], x0, opts);

    colorIdx = 1 + round((size(cmap,1)-1) * (cvals(k)-cmin) / max(cmax-cmin, eps));
    colorIdx = min(max(colorIdx,1), size(cmap,1));
    thisColor = cmap(colorIdx, :);

    plot(out.x(:,1), out.x(:,2), '-', 'Color', thisColor, 'LineWidth', 1.0);
end

colormap(cmap);
cb = colorbar;
ylabel(cb, cLabel, 'FontSize', 12, 'FontWeight', 'bold');
clim([cmin cmax]);

xlabel('x [-]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('y [-]', 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('JPL %s Family in %s CR3BP', ...
    capitalizeFirst(familyName), formatSystemName(sysName)), ...
    'FontSize', 15, 'FontWeight', 'bold');

% ----------------------------------------------------------
% Plot full family in 3D
% ----------------------------------------------------------
figure('Color','w');
hold on
grid on
box on
axis equal

plot3(-mu, 0, 0, 'bo', 'MarkerSize', 10, 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
plot3(1-mu, 0, 0, 'ko', 'MarkerSize', 7, 'LineWidth', 1.5, 'MarkerFaceColor', 'k');

for k = 1:nUse
    x0 = orbits(k).state;
    T = orbits(k).period;

    out = astro.propagators.propagate( ...
        @(t,x) astro.cr3bp.eomCR3BP(t, x, mu), ...
        [0, T], x0, opts);

    colorIdx = 1 + round((size(cmap,1)-1) * (cvals(k)-cmin) / max(cmax-cmin, eps));
    colorIdx = min(max(colorIdx,1), size(cmap,1));
    thisColor = cmap(colorIdx, :);

    plot3(out.x(:,1), out.x(:,2), out.x(:,3), '-', ...
        'Color', thisColor, 'LineWidth', 1.0);
end

colormap(cmap);
cb = colorbar;
ylabel(cb, cLabel, 'FontSize', 12, 'FontWeight', 'bold');
clim([cmin cmax]);

xlabel('x [-]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('y [-]', 'FontSize', 13, 'FontWeight', 'bold');
zlabel('z [-]', 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('JPL %s Family in %s CR3BP (3D)', ...
    capitalizeFirst(familyName), formatSystemName(sysName)), ...
    'FontSize', 15, 'FontWeight', 'bold');
view(35, 25)

% ----------------------------------------------------------
% Console diagnostics
% ----------------------------------------------------------
Cfirst = astro.cr3bp.jacobiConstant(orbits(1).state.', mu);
Clast  = astro.cr3bp.jacobiConstant(orbits(end).state.', mu);

fprintf('\nEndpoint family members:\n');
fprintf('  First orbit:  period = %.6f TU, Jacobi = %.6f, stability = %.6f\n', ...
    orbits(1).period, Cfirst, orbits(1).stability);
fprintf('  Last orbit:   period = %.6f TU, Jacobi = %.6f, stability = %.6f\n', ...
    orbits(end).period, Clast, orbits(end).stability);

fprintf('\nDone.\n');

function s = capitalizeFirst(str)
str = char(string(str));
if isempty(str)
    s = str;
else
    s = [upper(str(1)) str(2:end)];
end
end

function s = formatSystemName(sysName)
sysName = char(string(sysName));
s = strrep(sysName, '-', ' ');
parts = split(string(s));
for i = 1:numel(parts)
    p = char(parts(i));
    if ~isempty(p)
        parts(i) = string([upper(p(1)) p(2:end)]);
    end
end
s = char(join(parts, ' '));
end