function section = collectPoincareSection(seedStates, tMax, mu, xSection, direction)
%COLLECTPOINCARESECTION Propagate multiple seeds to a section x = xSection.
%
% INPUTS
%   seedStates : [6 x N] matrix of seed states
%   tMax       : max propagation time (sign gives forward/backward)
%   mu         : CR3BP mass parameter
%   xSection   : x-plane value
%   direction  : crossing direction for event
%
% OUTPUT
%   section : struct with fields
%       .points      [6 x M] crossing states
%       .times       [1 x M] crossing times
%       .hitMask     [1 x N]
%       .trajectories cell array of propagation outputs

[rows, N] = size(seedStates);
if rows ~= 6
    error('seedStates must be 6 x N.');
end

points = NaN(6, N);
times = NaN(1, N);
hitMask = false(1, N);
trajectories = cell(1, N);

for k = 1:N
    out = astro.cr3bp.propagateToSection(seedStates(:,k), tMax, mu, xSection, direction, true);
    trajectories{k} = out;

    if out.hit
        points(:,k) = out.xCross;
        times(k) = out.tCross;
        hitMask(k) = true;
    end
end

section = struct();
section.points = points(:, hitMask);
section.times = times(hitMask);
section.hitMask = hitMask;
section.trajectories = trajectories;
end