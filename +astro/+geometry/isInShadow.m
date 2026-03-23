function out = isInShadow(rSc, rSun, rBody, bodyRadius)
%ISINSHADOW Cylindrical shadow test for eclipse behind a spherical body.
%
% INPUTS
%   rSc        : spacecraft position [3x1 km]
%   rSun       : Sun position [3x1 km]
%   rBody      : occulting body center [3x1 km]
%   bodyRadius : occulting body radius [km]
%
% OUTPUT
%   out : struct with fields
%       .inShadow
%       .behindBody
%       .crossTrack
%       .alongTrack
%
% NOTES
%   This is a simple cylindrical shadow model:
%   - define the anti-Sun direction from the body
%   - project spacecraft relative position onto that axis
%   - if spacecraft is behind the body and within bodyRadius cross-track,
%     it is considered in shadow
%
%   This is a first-order eclipse model, useful for many preliminary studies.
%   Later this can be upgraded to umbra/penumbra conical geometry.

rSc = rSc(:);
rSun = rSun(:);
rBody = rBody(:);

sunDir = rSun - rBody;
sunDirNorm = norm(sunDir);

if sunDirNorm <= 0
    error('Sun and occulting body cannot have the same position.');
end

sHat = sunDir / sunDirNorm;          % toward Sun
aHat = -sHat;                        % anti-Sun shadow axis

rho = rSc - rBody;

alongTrack = dot(rho, aHat);         % positive if behind body
rhoPerp = rho - alongTrack * aHat;
crossTrack = norm(rhoPerp);

behindBody = alongTrack > 0;
inShadow = behindBody && (crossTrack <= bodyRadius);

out.inShadow = inShadow;
out.behindBody = behindBody;
out.crossTrack = crossTrack;
out.alongTrack = alongTrack;
end