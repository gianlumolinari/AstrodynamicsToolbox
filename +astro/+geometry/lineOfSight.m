function out = lineOfSight(r1, r2, rBody, bodyRadius)
%LINEOFSIGHT Checks whether the segment between two points is blocked by a spherical body.
%
% INPUTS
%   r1         : first point [3x1 km]
%   r2         : second point [3x1 km]
%   rBody      : body center [3x1 km]
%   bodyRadius : body radius [km]
%
% OUTPUT
%   out : struct with fields
%       .blocked         true if line segment intersects the body
%       .dMin            minimum distance from body center to infinite line [km]
%       .dMinSegment     minimum distance from body center to finite segment [km]
%       .closestPoint    closest point on the segment [3x1 km]
%       .tau             segment parameter in [0,1]
%
% NOTES
%   The line segment is parameterized as:
%       r(tau) = r1 + tau*(r2-r1),  0 <= tau <= 1
%
%   If the closest point on the segment lies within the body radius, then
%   the line of sight is blocked.

r1 = r1(:);
r2 = r2(:);
rBody = rBody(:);

dr = r2 - r1;
L2 = dot(dr, dr);

if L2 <= 0
    error('lineOfSight requires r1 and r2 to be distinct points.');
end

tauLine = dot(rBody - r1, dr) / L2;
tau = min(max(tauLine, 0), 1);

closestPoint = r1 + tau * dr;

dMin = norm(cross(dr, rBody - r1)) / norm(dr);
dMinSegment = norm(closestPoint - rBody);

blocked = (dMinSegment <= bodyRadius);

out.blocked = blocked;
out.dMin = dMin;
out.dMinSegment = dMinSegment;
out.closestPoint = closestPoint;
out.tau = tau;
end