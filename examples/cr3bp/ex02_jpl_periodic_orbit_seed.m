clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/AstroToolbox')
startup

% ==========================================================
% JPL periodic-orbit seed example
%
% Queries a JPL CR3BP family, parses one orbit, propagates it
% with our own CR3BP equations, and compares Jacobi constant.
% ==========================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('JPL Periodic Orbit Seed Demo\n');
fprintf('============================================================\n');

% ----------------------------------------------------------
% Query JPL for an Earth-Moon northern L1 halo family subset
% ----------------------------------------------------------
data = astro.cr3bp.queryJPLPeriodicOrbits( ...
    'sys', 'earth-moon', ...
    'family', 'halo', ...
    'libr', 1, ...
    'branch', 'N', ...
    'periodmax', 2.2, ...
    'periodunits', 'TU');

fprintf('\nJPL API returned:\n');
fprintf('  System     : %s\n', data.system.name);
fprintf('  Family     : %s\n', data.family);
if isfield(data, 'libration_point')
    fprintf('  Libration  : %s\n', data.libration_point);
end
if isfield(data, 'branch')
    fprintf('  Branch     : %s\n', data.branch);
end
fprintf('  Count      : %s\n', data.count);
fprintf('  API ver.   : %s\n', data.signature.version);

% ----------------------------------------------------------
% Pick one orbit from the returned family
% ----------------------------------------------------------
idx = 1;
orbit = astro.cr3bp.parseJPLPeriodicOrbit(data, idx);

fprintf('\nSelected orbit row %d:\n', idx);
fprintf('  x,y,z       = [%.12f, %.12f, %.12f]\n', orbit.state(1), orbit.state(2), orbit.state(3));
fprintf('  vx,vy,vz    = [%.12f, %.12f, %.12f]\n', orbit.state(4), orbit.state(5), orbit.state(6));
fprintf('  Jacobi      : %.12f\n', orbit.jacobi);
fprintf('  Period (TU) : %.12f\n', orbit.period);
fprintf('  Stability   : %.12f\n', orbit.stability);
fprintf('  mu          : %.15f\n', orbit.mu);
fprintf('  lunit [km]  : %.6f\n', orbit.lunit_km);
fprintf('  tunit [s]   : %.6f\n', orbit.tunit_s);

% ----------------------------------------------------------
% Propagate the seed orbit in our CR3BP model
% ----------------------------------------------------------
x0 = orbit.state;
mu = orbit.mu;

opts.RelTol = 1e-12;
opts.AbsTol = 1e-12;
opts.Solver = 'ode113';

out = astro.propagators.propagate( ...
    @(t,x) astro.cr3bp.eomCR3BP(t, x, mu), ...
    [0, orbit.period], x0, opts);

C0 = astro.cr3bp.jacobiConstant(x0.', mu);
Cend = astro.cr3bp.jacobiConstant(out.x(end,:), mu);

fprintf('\nCR3BP propagation diagnostics:\n');
fprintf('  Initial Jacobi       : %.12f\n', C0);
fprintf('  Final Jacobi         : %.12f\n', Cend);
fprintf('  Jacobi drift         : %.6e\n', abs(Cend - C0));

stateErr = norm(out.x(end,:).' - x0);
fprintf('  Periodicity defect   : %.6e\n', stateErr);

% ----------------------------------------------------------
% Plot x-y projection
% ----------------------------------------------------------
figure('Color','w');
hold on
axis equal
grid on
box on

plot(out.x(:,1), out.x(:,2), 'b-', 'LineWidth', 1.4);

% Primaries
plot(-mu, 0, 'bo', 'MarkerSize', 10, 'LineWidth', 1.4, 'MarkerFaceColor', 'b');
plot(1-mu, 0, 'ko', 'MarkerSize', 7, 'LineWidth', 1.4, 'MarkerFaceColor', 'k');

% Lagrange points
L = astro.cr3bp.lagrangePoints(mu);
plot(L.L1(1), L.L1(2), 'rs', 'MarkerSize', 7, 'LineWidth', 1.4);
plot(L.L2(1), L.L2(2), 'rs', 'MarkerSize', 7, 'LineWidth', 1.4);

text(L.L1(1), L.L1(2), '  L1', 'FontSize', 10, 'FontWeight', 'bold');
text(L.L2(1), L.L2(2), '  L2', 'FontSize', 10, 'FontWeight', 'bold');

xlabel('x [-]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('y [-]', 'FontSize', 13, 'FontWeight', 'bold');
title('JPL Halo-Orbit Seed Propagated in Our CR3BP Model', ...
    'FontSize', 15, 'FontWeight', 'bold');
legend('Periodic orbit', 'Primary', 'Secondary', 'L1/L2', 'Location', 'best');

% ----------------------------------------------------------
% 3D plot
% ----------------------------------------------------------
figure('Color','w');
hold on
grid on
box on
axis equal

plot3(out.x(:,1), out.x(:,2), out.x(:,3), 'b-', 'LineWidth', 1.4);
plot3(-mu, 0, 0, 'bo', 'MarkerSize', 10, 'LineWidth', 1.4, 'MarkerFaceColor', 'b');
plot3(1-mu, 0, 0, 'ko', 'MarkerSize', 7, 'LineWidth', 1.4, 'MarkerFaceColor', 'k');

xlabel('x [-]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('y [-]', 'FontSize', 13, 'FontWeight', 'bold');
zlabel('z [-]', 'FontSize', 13, 'FontWeight', 'bold');
title('JPL Halo-Orbit Seed: 3D View', 'FontSize', 15, 'FontWeight', 'bold');
view(35, 25)

fprintf('\nDone.\n');