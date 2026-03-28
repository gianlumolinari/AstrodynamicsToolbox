function h = plotManifoldFamily3D(atlas, opts)
%PLOTMANIFOLDFAMILY3D Plot a 3D manifold family.
%
% INPUTS
%   atlas : manifold atlas
%   opts  : struct with optional fields
%       .ax
%       .Color
%       .LineWidth
%       .ShowOrbit
%       .OrbitColor
%       .OrbitLineWidth
%       .ShowPrimaries
%       .mu
%       .ShowLibrationPoints
%       .View
%
% OUTPUT
%   h : struct with plot handles

if nargin < 2 || isempty(opts)
    opts = struct();
end

if ~isfield(opts,'ax'), opts.ax = []; end
if ~isfield(opts,'Color'), opts.Color = [0 0 1]; end
if ~isfield(opts,'LineWidth'), opts.LineWidth = 0.8; end
if ~isfield(opts,'ShowOrbit'), opts.ShowOrbit = true; end
if ~isfield(opts,'OrbitColor'), opts.OrbitColor = [0 0 0]; end
if ~isfield(opts,'OrbitLineWidth'), opts.OrbitLineWidth = 2.0; end
if ~isfield(opts,'ShowPrimaries'), opts.ShowPrimaries = true; end
if ~isfield(opts,'mu'), opts.mu = atlas.mu; end
if ~isfield(opts,'ShowLibrationPoints'), opts.ShowLibrationPoints = false; end
if ~isfield(opts,'View'), opts.View = [35 25]; end

if isempty(opts.ax)
    figure;
    ax = axes;
else
    ax = opts.ax;
end

hold(ax,'on');
grid(ax,'on');
box(ax,'on');

branches = atlas.branches;
h.branch = gobjects(numel(branches),1);

for k = 1:numel(branches)
    X = branches(k).x;
    h.branch(k) = plot3(ax, X(:,1), X(:,2), X(:,3), ...
        'Color', opts.Color, ...
        'LineWidth', opts.LineWidth);
end

h.orbit = gobjects(1);
if opts.ShowOrbit
    Xo = atlas.orbit.x;
    h.orbit = plot3(ax, Xo(:,1), Xo(:,2), Xo(:,3), ...
        'Color', opts.OrbitColor, ...
        'LineWidth', opts.OrbitLineWidth);
end

h.primaries = gobjects(0);
if opts.ShowPrimaries
    mu = opts.mu;
    h.primaries(1) = plot3(ax, -mu, 0, 0, 'ko', ...
        'MarkerFaceColor','k', 'MarkerSize', 6);
    h.primaries(2) = plot3(ax, 1-mu, 0, 0, 'o', ...
        'Color',[0.4 0.4 0.4], ...
        'MarkerFaceColor',[0.7 0.7 0.7], ...
        'MarkerSize', 8);
end

h.libration = gobjects(0);
if opts.ShowLibrationPoints
    L = astro.cr3bp.lagrangePoints(opts.mu);
    h.libration(1) = plot3(ax, L.L1(1), L.L1(2), L.L1(3), 'yo', ...
        'MarkerFaceColor','y', 'MarkerSize', 7);
    h.libration(2) = plot3(ax, L.L2(1), L.L2(2), L.L2(3), 'co', ...
        'MarkerFaceColor','c', 'MarkerSize', 7);
end

axis(ax,'equal');
xlabel(ax,'x');
ylabel(ax,'y');
zlabel(ax,'z');
view(ax, opts.View);
end