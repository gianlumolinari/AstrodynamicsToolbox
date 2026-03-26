function h = plotManifoldFamily2D(atlas, opts)
%PLOTMANIFOLDFAMILY2D Plot a 2D manifold family with custom styling.
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

    if nargin < 2 || isempty(opts)
        opts = struct();
    end

    if ~isfield(opts,'ax'), opts.ax = []; end
    if ~isfield(opts,'Color'), opts.Color = [0 0 1]; end
    if ~isfield(opts,'LineWidth'), opts.LineWidth = 0.8; end
    if ~isfield(opts,'ShowOrbit'), opts.ShowOrbit = false; end
    if ~isfield(opts,'OrbitColor'), opts.OrbitColor = [0 0 0]; end
    if ~isfield(opts,'OrbitLineWidth'), opts.OrbitLineWidth = 2.0; end
    if ~isfield(opts,'ShowPrimaries'), opts.ShowPrimaries = false; end

    if isempty(opts.ax)
        figure;
        ax = axes;
    else
        ax = opts.ax;
    end

    hold(ax,'on');
    grid(ax,'on');

    branches = atlas.branches;
    h.branch = gobjects(numel(branches),1);

    for k = 1:numel(branches)
        X = branches(k).x;
        h.branch(k) = plot(ax, X(:,1), X(:,2), ...
            'Color', opts.Color, ...
            'LineWidth', opts.LineWidth);
    end

    h.orbit = gobjects(1);
    if opts.ShowOrbit
        Xo = atlas.orbit.x;
        h.orbit = plot(ax, Xo(:,1), Xo(:,2), ...
            'Color', opts.OrbitColor, ...
            'LineWidth', opts.OrbitLineWidth);
    end

    h.primaries = gobjects(0);
    if opts.ShowPrimaries
        mu = atlas.mu;
        h.primaries(1) = plot(ax, -mu, 0, 'ko', ...
            'MarkerFaceColor','k', 'MarkerSize', 6);
        h.primaries(2) = plot(ax, 1-mu, 0, 'o', ...
            'Color',[0.4 0.4 0.4], ...
            'MarkerFaceColor',[0.7 0.7 0.7], ...
            'MarkerSize', 8);
    end

    axis(ax,'equal');
    xlabel(ax,'x');
    ylabel(ax,'y');
end