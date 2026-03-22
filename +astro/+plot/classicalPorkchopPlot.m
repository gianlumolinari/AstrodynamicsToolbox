function ax = classicalPorkchopPlot(depDates, arrDates, metricGrid, metricLevels, ...
    tofGrid, tofLevels, plotTitle, metricLabel, bestPoint)
%CLASSICALPORKCHOPPLOT Creates a classical contour-style porkchop plot.

if nargin < 9
    bestPoint = [];
end

depNum = datenum(depDates);
arrNum = datenum(arrDates);

[DEP, ARR] = meshgrid(depNum, arrNum);

figureColor = 'w';
axesColor = 'w';

set(gcf, 'Color', figureColor)
set(gca, 'Color', axesColor)

hold on
box on
grid on

% Colored contours for main metric
[Cm, hm] = contour(DEP, ARR, metricGrid, metricLevels, 'LineWidth', 1.3);
clabel(Cm, hm, 'FontSize', 8, 'Color', 'k')

colormap(turbo)
caxis([min(metricLevels) max(metricLevels)])
cb = colorbar;
ylabel(cb, metricLabel)

% Black TOF contours
if ~isempty(tofGrid)
    [Ct, ht] = contour(DEP, ARR, tofGrid, tofLevels, 'k', 'LineWidth', 1.5);
    clabel(Ct, ht, 'FontSize', 8, 'Color', 'k', 'FontWeight', 'bold')
end

% Best point marker
if ~isempty(bestPoint)
    xBest = datenum(bestPoint.depDate);
    yBest = datenum(bestPoint.arrDate);
    plot(xBest, yBest, 'ko', 'MarkerFaceColor', 'w', ...
        'MarkerSize', 7, 'LineWidth', 1.5)
    
    if isfield(bestPoint, 'label') && ~isempty(bestPoint.label)
        text(xBest, yBest, ['  ' bestPoint.label], ...
            'Color', 'k', 'FontWeight', 'bold', 'FontSize', 9)
    end
end

xlabel('Departure Date')
ylabel('Arrival Date')
title(plotTitle)

datetick('x', 'dd-mmm-yyyy', 'keeplimits')
datetick('y', 'dd-mmm-yyyy', 'keeplimits')
xtickangle(30)

ax = gca;
end