function out = isVisibleFromStation(rSc, rStation, minElevationDeg)
%ISVISIBLEFROMSTATION Visibility test above an elevation mask.
%
% INPUTS
%   rSc             : spacecraft position [3x1 km]
%   rStation        : station position [3x1 km]
%   minElevationDeg : minimum elevation mask [deg]
%
% OUTPUT
%   out : struct with fields
%       .isVisible
%       .elevation
%
% NOTE
%   This is a pure geometric visibility test.

if nargin < 3 || isempty(minElevationDeg)
    minElevationDeg = 0;
end

elevation = astro.visibility.elevationAngle(rSc, rStation);
isVisible = elevation.elevationDeg >= minElevationDeg;

out.isVisible = isVisible;
out.elevation = elevation;
end