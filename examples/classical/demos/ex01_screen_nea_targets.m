clc;
clear;
close all;

cd('/Users/gianlucamolinari/Desktop/AstroToolbox')
startup

% ==========================================================
% Small-body screening demo
%
% Purpose:
%   quick triage of a few well-known NEA targets using JPL SBDB
% ==========================================================

targets = {
    'Apophis'
    'Bennu'
    'Eros'
    'Itokawa'
    'Ryugu'
    'Didymos'
};

fprintf('\n');
fprintf('============================================================\n');
fprintf('NEA Screening Demo (JPL SBDB)\n');
fprintf('============================================================\n');

results = astro.smallbody.screenNEOTargets(targets);

fprintf('\nQueried targets:\n');
for k = 1:numel(results)
    fprintf('  %-10s  a = %8.4f au   e = %7.4f   i = %7.3f deg   MOID = %8.5f au   H = %6.2f\n', ...
        results(k).query, ...
        results(k).a_au, ...
        results(k).e, ...
        results(k).i_deg, ...
        results(k).moid_au, ...
        results(k).H);
end

% ----------------------------------------------------------
% Simple mission-relevance score
% Lower is better
% ----------------------------------------------------------
score = NaN(numel(results),1);
for k = 1:numel(results)
    if ~isnan(results(k).a_au) && ~isnan(results(k).e) && ...
       ~isnan(results(k).i_deg) && ~isnan(results(k).moid_au)

        score(k) = ...
            2.0 * abs(results(k).a_au - 1.0) + ...
            1.5 * results(k).e + ...
            0.03 * results(k).i_deg + ...
            20.0 * results(k).moid_au;
    end
end

[scoreSorted, I] = sort(score, 'ascend', 'MissingPlacement', 'last');
resultsSorted = results(I);

fprintf('\nRanked shortlist (heuristic score; lower = more accessible-like):\n');
for k = 1:numel(resultsSorted)
    fprintf('  %d) %-10s  score = %8.4f\n', ...
        k, resultsSorted(k).query, scoreSorted(k));
end

% ----------------------------------------------------------
% Plot: a vs e
% ----------------------------------------------------------
figure('Color','w');
hold on
grid on
box on

aVals = [results.a_au];
eVals = [results.e];

plot(aVals, eVals, 'o', 'MarkerSize', 8, 'LineWidth', 1.3);

for k = 1:numel(results)
    if ~isnan(results(k).a_au) && ~isnan(results(k).e)
        text(results(k).a_au, results(k).e, ['  ' results(k).query], ...
            'FontSize', 10, 'FontWeight', 'bold');
    end
end

xlabel('Semi-major axis a [au]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Eccentricity e [-]', 'FontSize', 13, 'FontWeight', 'bold');
title('Candidate NEA Screening: a-e Map', 'FontSize', 15, 'FontWeight', 'bold');

fprintf('\nDone.\n');