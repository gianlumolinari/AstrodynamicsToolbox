function ax = classicalPorkchopPlot(depDates, arrDates, metricGrid, metricLevels, ...
    tofGrid, tofLevels, plotTitle, metricLabel, bestPoint)
%CLASSICALPORKCHOPPLOT Smoothed classical contour-style porkchop plot.
%
% INPUTS
%   depDates     : vector of departure datetimes
%   arrDates     : vector of arrival datetimes
%   metricGrid   : matrix [numArr x numDep]
%   metricLevels : contour levels for metric
%   tofGrid      : matrix [numArr x numDep], time of flight in days
%   tofLevels    : contour levels for TOF
%   plotTitle    : title string
%   metricLabel  : colorbar label
%   bestPoint    : optional struct with fields:
%                  .depDate
%                  .arrDate
%                  .label
%
% OUTPUT
%   ax : axes handle

if nargin < 9
    bestPoint = [];
end

depNum = datenum(depDates);
arrNum = datenum(arrDates);

[DEP, ARR] = meshgrid(depNum, arrNum);

% ----------------------------------------------------------
% Interpolate to a finer mesh for smoother classical contours
% ----------------------------------------------------------
depFine = linspace(min(depNum), max(depNum), 4*numel(depNum));
arrFine = linspace(min(arrNum), max(arrNum), 4*numel(arrNum));
[DEPF, ARRF] = meshgrid(depFine, arrFine);

metricFine = interp2(DEP, ARR, metricGrid, DEPF, ARRF, 'linear');
tofFine    = interp2(DEP, ARR, tofGrid,    DEPF, ARRF, 'linear');

% ----------------------------------------------------------
% Plot
% ----------------------------------------------------------
set(gcf, 'Color', 'w')
set(gca, 'Color', 'w')
hold on
box on
grid on

% Light filled background
contourf(DEPF, ARRF, metricFine, 100, 'LineStyle', 'none');
colormap(turbo)
cb = colorbar;
ylabel(cb, metricLabel)

% Colored metric contours
[Cm, hm] = contour(DEPF, ARRF, metricFine, metricLevels, 'LineWidth', 1.2);
clabel(Cm, hm, 'FontSize', 8, 'Color', 'k')

% Black TOF contours
[Ct, ht] = contour(DEPF, ARRF, tofFine, tofLevels, 'k', 'LineWidth', 1.4);
clabel(Ct, ht, 'FontSize', 8, 'Color', 'k', 'FontWeight', 'bold')

% Best point
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