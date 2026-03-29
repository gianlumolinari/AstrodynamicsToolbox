clc;
clear;
close all;

% INSPECT_VALIDATED_ORBIT_DATABASE
% Browse the current validated CR3BP orbit database and visualize saved orbits.

systemName = 'earth_moon';
mu = astro.cr3bp.getSystemMu(systemName);

files = astro.periodic.listValidatedOrbits(systemName);

if isempty(files)
    error('No validated orbit files found.');
end

fprintf('\nLoading validated orbit database...\n');

db = struct([]);
for k = 1:numel(files)
    S = load(fullfile(files(k).folder, files(k).name));
    S.fileName = files(k).name;
    db(k,1) = S;
end

fprintf('\nValidated orbit summary\n');
fprintf('-----------------------\n');
for k = 1:numel(db)
    fprintf('%2d) %-35s  T = %.10f   closure = %.3e\n', ...
        k, db(k).fileName, db(k).T, db(k).closureError);
end

fig = figure;
ax = axes(fig);
hold(ax,'on');
grid(ax,'on');
box(ax,'on');

legendHandles = gobjects(numel(db),1);
legendLabels = cell(numel(db),1);

for k = 1:numel(db)
    orbit = astro.periodic.buildOrbitStruct(db(k).state0, db(k).T, mu, ...
        struct('family', db(k).family, ...
               'libr', db(k).libr, ...
               'system', systemName, ...
               'source', db(k).source, ...
               'dimension', localInferDimension(db(k).family)), ...
        struct('verbose', false));

    if strcmpi(localInferDimension(db(k).family), '2D')
        legendHandles(k) = plot3(ax, orbit.x(:,1), orbit.x(:,2), zeros(size(orbit.x,1),1), ...
            'LineWidth', 1.8);
    else
        legendHandles(k) = plot3(ax, orbit.x(:,1), orbit.x(:,2), orbit.x(:,3), ...
            'LineWidth', 1.8);
    end
    legendLabels{k} = db(k).fileName;
end

L = astro.cr3bp.lagrangePoints(mu);
plot3(ax, -mu, 0, 0, 'ko', 'MarkerFaceColor','k', 'MarkerSize', 6);
plot3(ax, 1-mu, 0, 0, 'o', ...
    'Color',[0.4 0.4 0.4], ...
    'MarkerFaceColor',[0.7 0.7 0.7], ...
    'MarkerSize', 8);
plot3(ax, L.L1(1), L.L1(2), L.L1(3), 'yo', 'MarkerFaceColor','y', 'MarkerSize', 7);
plot3(ax, L.L2(1), L.L2(2), L.L2(3), 'co', 'MarkerFaceColor','c', 'MarkerSize', 7);

axis(ax,'equal');
xlabel(ax,'x');
ylabel(ax,'y');
zlabel(ax,'z');
title(ax,'Validated CR3BP orbit database');
view(ax,[35 25]);

legend(ax, legendHandles, legendLabels, 'Interpreter','none', 'Location','bestoutside');

function dim = localInferDimension(familyName)
    switch lower(string(familyName))
        case "planar_lyapunov"
            dim = '2D';
        otherwise
            dim = '3D';
    end
end