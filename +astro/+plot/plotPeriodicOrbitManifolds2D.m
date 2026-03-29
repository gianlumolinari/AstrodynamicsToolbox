function h = plotPeriodicOrbitManifolds2D(man, orbit, mu, plotTitle, ax)
%PLOTPERIODICORBITMANIFOLDS2D Plot stable/unstable manifolds for a planar periodic orbit.
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

h.unstable = astro.plot.plotManifoldFamily2D(man.unstable, struct( ...
    'ax', ax, ...
    'Color', [1 0 0], ...
    'LineWidth', 0.8, ...
    'ShowOrbit', false, ...
    'ShowPrimaries', false));

h.stable = astro.plot.plotManifoldFamily2D(man.stable, struct( ...
    'ax', ax, ...
    'Color', [0 0 1], ...
    'LineWidth', 0.8, ...
    'ShowOrbit', false, ...
    'ShowPrimaries', false));

h.orbit = plot(ax, orbit.x(:,1), orbit.x(:,2), 'k', 'LineWidth', 2.0);

L = astro.cr3bp.lagrangePoints(mu);
h.earth = plot(ax, -mu, 0, 'ko', 'MarkerFaceColor','k', 'MarkerSize', 6);
h.moon  = plot(ax, 1-mu, 0, 'o', ...
    'Color',[0.4 0.4 0.4], ...
    'MarkerFaceColor',[0.7 0.7 0.7], ...
    'MarkerSize', 8);
h.L1 = plot(ax, L.L1(1), L.L1(2), 'yo', 'MarkerFaceColor','y', 'MarkerSize', 7);
h.L2 = plot(ax, L.L2(1), L.L2(2), 'co', 'MarkerFaceColor','c', 'MarkerSize', 7);

axis(ax,'equal');
xlabel(ax,'x');
ylabel(ax,'y');
title(ax, plotTitle, 'Interpreter','none');

legend(ax, ...
    [h.unstable.branch(1), h.stable.branch(1), h.orbit, h.earth, h.moon, h.L1, h.L2], ...
    {'Unstable manifold','Stable manifold','Periodic orbit','Earth','Moon','L_1','L_2'}, ...
    'Location','best');
end