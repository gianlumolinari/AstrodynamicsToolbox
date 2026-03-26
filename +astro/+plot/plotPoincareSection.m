function h = plotPoincareSection(sectionData, opts)
%PLOTPOINCARESECTION Plot Poincare section points.
%
% INPUTS
%   sectionData : struct with fields .x and optionally .label
%   opts        : plotting options
%
% OUTPUT
%   h : graphics handle

    arguments
        sectionData struct
        opts.ax = []
        opts.xIndex (1,1) double {mustBeInteger,mustBePositive} = 2
        opts.yIndex (1,1) double {mustBeInteger,mustBePositive} = 5
    end

    if isempty(opts.ax)
        figure;
        ax = axes;
    else
        ax = opts.ax;
    end

    hold(ax,'on');
    grid(ax,'on');

    h = plot(ax, sectionData.x(:,opts.xIndex), sectionData.x(:,opts.yIndex), '.');

    xlabel(ax, sprintf('x_%d', opts.xIndex));
    ylabel(ax, sprintf('x_%d', opts.yIndex));
    title(ax, 'Poincare section');
end