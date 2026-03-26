function h = plotManifoldTube(atlas, opts)
%PLOTMANIFOLDTUBE Plot a manifold atlas as a set of strands.

    if nargin < 2 || isempty(opts)
        opts = struct();
    end

    if ~isfield(opts,'ax'), opts.ax = []; end
    if ~isfield(opts,'dimension'), opts.dimension = '2d'; end
    if ~isfield(opts,'showOrbit'), opts.showOrbit = true; end
    if ~isfield(opts,'showPrimaries'), opts.showPrimaries = true; end
    if ~isfield(opts,'lineWidth'), opts.lineWidth = 0.8; end
    if ~isfield(opts,'orbitLineWidth'), opts.orbitLineWidth = 2.0; end

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

        if strcmpi(opts.dimension,'2d')
            h.branch(k) = plot(ax, X(:,1), X(:,2), 'LineWidth', opts.lineWidth);
        else
            h.branch(k) = plot3(ax, X(:,1), X(:,2), X(:,3), 'LineWidth', opts.lineWidth);
        end
    end

    h.orbit = gobjects(1);
    if opts.showOrbit
        Xo = atlas.orbit.x;
        if strcmpi(opts.dimension,'2d')
            h.orbit = plot(ax, Xo(:,1), Xo(:,2), 'k', 'LineWidth', opts.orbitLineWidth);
        else
            h.orbit = plot3(ax, Xo(:,1), Xo(:,2), Xo(:,3), 'k', 'LineWidth', opts.orbitLineWidth);
        end
    end

    h.primaries = gobjects(0);
    if opts.showPrimaries
        mu = atlas.mu;

        if strcmpi(opts.dimension,'2d')
            h.primaries(1) = plot(ax, -mu, 0, 'o', 'MarkerSize', 8, 'LineWidth', 1.5);
            h.primaries(2) = plot(ax, 1-mu, 0, 'o', 'MarkerSize', 6, 'LineWidth', 1.5);
        else
            h.primaries(1) = plot3(ax, -mu, 0, 0, 'o', 'MarkerSize', 8, 'LineWidth', 1.5);
            h.primaries(2) = plot3(ax, 1-mu, 0, 0, 'o', 'MarkerSize', 6, 'LineWidth', 1.5);
        end
    end

    axis(ax,'equal');
    xlabel(ax,'x');
    ylabel(ax,'y');
    if strcmpi(opts.dimension,'3d')
        zlabel(ax,'z');
        view(ax,3);
    end

    title(ax, sprintf('%s manifold tube', atlas.type));
end