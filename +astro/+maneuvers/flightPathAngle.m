function gamma = flightPathAngle(e, theta)
%FLIGHTPATHANGLE Computes Keplerian flight-path angle.
%
% INPUTS
%   e     : eccentricity
%   theta : true anomaly [rad]
%
% OUTPUT
%   gamma : flight-path angle [rad]
%
% RELATION
%   tan(gamma) = e*sin(theta) / (1 + e*cos(theta))

gamma = atan2(e .* sin(theta), 1 + e .* cos(theta));
end