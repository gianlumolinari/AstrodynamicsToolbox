function out = isOcculted(rObserver, rTarget, rOccultingBody, occultingRadius)
%ISOCCULTED Checks whether a target is occulted from an observer by a spherical body.
%
% INPUTS
%   rObserver        : observer position [3x1 km]
%   rTarget          : target position [3x1 km]
%   rOccultingBody   : occulting body center [3x1 km]
%   occultingRadius  : occulting body radius [km]
%
% OUTPUT
%   out : struct with fields
%       .isOcculted
%       .los
%
% NOTES
%   This is a wrapper around astro.geometry.lineOfSight.

los = astro.geometry.lineOfSight(rObserver, rTarget, rOccultingBody, occultingRadius);

out.isOcculted = los.blocked;
out.los = los;
end