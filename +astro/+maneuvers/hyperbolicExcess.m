function vInfVec = hyperbolicExcess(vScHelio, vPlanetHelio)
%HYPERBOLICEXCESS Computes hyperbolic excess velocity vector.
%
% INPUTS
%   vScHelio     : spacecraft heliocentric velocity [3x1 km/s]
%   vPlanetHelio : planet heliocentric velocity [3x1 km/s]
%
% OUTPUT
%   vInfVec      : hyperbolic excess velocity vector wrt planet [3x1 km/s]

vInfVec = vScHelio(:) - vPlanetHelio(:);
end