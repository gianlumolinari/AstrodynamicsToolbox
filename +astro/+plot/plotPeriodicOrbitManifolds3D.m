function h = plotPeriodicOrbitManifolds3D(man, orbit, mu, plotTitle, ax)
%PLOTPERIODICORBITMANIFOLDS3D Plot stable/unstable manifolds for a 3D periodic orbit.
%
% INPUTS
%   man       : output of astro.manifolds.generatePeriodicOrbitManifolds
%   orbit     : periodic orbit struct
%   mu        : CR3BP mass parameter
%   plotTitle : figure/axes title
%   ax        : optional axes handle
%
% OUTPUT
%   h : struct with graphics handles

if nargin < 4 || isempty(plotTitle)
    plotTitle = 'Periodic orbit manifolds';
end

if nargin < 5 || isempty(ax)
    fig = figure;
    ax = axes(fig);
end

hold(ax,'on');
grid(ax,'on');
box(ax,'on');

h = struct();

h.unstable = astro.plot.plotManifoldFamily3D(man.unstable, struct( ...
    'ax', ax, ...
    'Color', [1 0 0], ...
    'LineWidth', 0.8, ...
    'ShowOrbit', false, ...
    'ShowPrimaries', false, ...
    'ShowLibrationPoints', false));

h.stable = astro.plot.plotManifoldFamily3D(man.stable, struct( ...
    'ax', ax, ...
    'Color', [0 0 1], ...
    'LineWidth', 0.8, ...
    'ShowOrbit', false, ...
    'ShowPrimaries', false, ...
    'ShowLibrationPoints', false));

h.orbit = plot3(ax, orbit.x(:,1), orbit.x(:,2), orbit.x(:,3), 'k', 'LineWidth', 2.0);

L = astro.cr3bp.lagrangePoints(mu);
h.earth = plot3(ax, -mu, 0, 0, 'ko', 'MarkerFaceColor','k', 'MarkerSize', 6);
h.moon  = plot3(ax, 1-mu, 0, 0, 'o', ...
    'Color',[0.4 0.4 0.4], ...
    'MarkerFaceColor',[0.7 0.7 0.7], ...
    'MarkerSize', 8);
h.L1 = plot3(ax, L.L1(1), L.L1(2), L.L1(3), 'yo', ...
    'MarkerFaceColor','y', 'MarkerSize', 7);
h.L2 = plot3(ax, L.L2(1), L.L2(2), L.L2(3), 'co', ...
    'MarkerFaceColor','c', 'MarkerSize', 7);

axis(ax,'equal');
xlabel(ax,'x');
ylabel(ax,'y');
zlabel(ax,'z');
title(ax, plotTitle, 'Interpreter','none');
view(ax,[35 25]);

legend(ax, ...
    [h.unstable.branch(1), h.stable.branch(1), h.orbit, h.earth, h.moon, h.L1, h.L2], ...
    {'Unstable manifold','Stable manifold','Periodic orbit','Earth','Moon','L_1','L_2'}, ...
    'Location','best');
end