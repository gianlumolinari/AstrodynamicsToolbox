function h = plotPlanetOrbit2D(a, varargin)
%PLOTPLANETORBIT2D Plots a circular 2D orbit of radius a centered at origin.
%
% INPUT
%   a : orbit radius [km]
%
% OPTIONAL INPUT
%   Additional plotting options
%
% OUTPUT
%   h : line handle

theta = linspace(0, 2*pi, 500);
x = a * cos(theta);
y = a * sin(theta);

h = plot(x, y, varargin{:});
grid on
axis equal
xlabel('x [km]')
ylabel('y [km]')
end