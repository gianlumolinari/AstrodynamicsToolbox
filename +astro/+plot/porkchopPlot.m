function ax = porkchopPlot(depDates, arrDates, Z, plotTitle, colorLabel, contourData, contourLevels, bestPoint)
%PORKCHOPPLOT Creates a colored porkchop plot with optional contours and marker.
%
% INPUTS
%   depDates      : vector of departure datetimes
%   arrDates      : vector of arrival datetimes
%   Z             : matrix of metric values [numArr x numDep]
%   plotTitle     : plot title
%   colorLabel    : colorbar label
%   contourData   : optional matrix for contour overlay [numArr x numDep]
%   contourLevels : optional contour levels
%   bestPoint     : optional struct with fields:
%                   .depDate
%                   .arrDate
%                   .label
%
% OUTPUT
%   ax : axes handle

if nargin < 6
    contourData = [];
end
if nargin < 7
    contourLevels = [];
end
if nargin < 8
    bestPoint = [];
end

depNum = datenum(depDates);
arrNum = datenum(arrDates);

imagesc(depNum, arrNum, Z);
set(gca, 'YDir', 'normal')
hold on
grid on

xlabel('Departure date')
ylabel('Arrival date')
title(plotTitle)

cb = colorbar;
ylabel(cb, colorLabel)

datetick('x', 'dd-mmm-yyyy', 'keeplimits')
datetick('y', 'dd-mmm-yyyy', 'keeplimits')
xtickangle(30)

% Optional contour overlay
if ~isempty(contourData)
    if isempty(contourLevels)
        [C, hC] = contour(depNum, arrNum, contourData, 'k');
    else
        [C, hC] = contour(depNum, arrNum, contourData, contourLevels, 'k');
    end
    clabel(C, hC, 'Color', 'k', 'FontSize', 8)
end

% Optional best-point marker
if ~isempty(bestPoint)
    xBest = datenum(bestPoint.depDate);
    yBest = datenum(bestPoint.arrDate);
    plot(xBest, yBest, 'wo', 'MarkerSize', 7, 'LineWidth', 1.5)
    
    if isfield(bestPoint, 'label') && ~isempty(bestPoint.label)
        text(xBest, yBest, ['  ' bestPoint.label], ...
            'Color', 'w', 'FontWeight', 'bold', 'FontSize', 9)
    end
end

ax = gca;
end