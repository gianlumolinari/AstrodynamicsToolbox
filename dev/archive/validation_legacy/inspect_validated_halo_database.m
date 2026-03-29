clc;
clear;
close all;

files = astro.periodic.loadValidatedFamily('earth_moon', 'halo', 1, 'north');

mu = 0.012150585609624;

fig = figure;
ax = axes(fig);
hold(ax,'on');
grid(ax,'on');
box(ax,'on');

for k = 1:numel(files)
    orbit = astro.periodic.buildOrbitStruct(files(k).state0, files(k).T, mu, ...
        struct('family', files(k).family, ...
               'libr', files(k).libr, ...
               'system', 'earth-moon', ...
               'source', files(k).source, ...
               'dimension', '3D'), ...
        struct('verbose', false));

    plot3(ax, orbit.x(:,1), orbit.x(:,2), orbit.x(:,3), 'LineWidth', 1.5);
end

L = astro.cr3bp.lagrangePoints(mu);
plot3(ax, -mu, 0, 0, 'ko', 'MarkerFaceColor','k', 'MarkerSize', 6);
plot3(ax, 1-mu, 0, 0, 'o', ...
    'Color',[0.4 0.4 0.4], ...
    'MarkerFaceColor',[0.7 0.7 0.7], ...
    'MarkerSize', 8);
plot3(ax, L.L1(1), L.L1(2), L.L1(3), 'yo', 'MarkerFaceColor','y', 'MarkerSize', 7);

axis(ax,'equal');
xlabel(ax,'x');
ylabel(ax,'y');
zlabel(ax,'z');
title(ax,'Validated halo L1 north database');
view(ax,[35 25]);