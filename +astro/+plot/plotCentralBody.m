function h = plotCentralBody(body)
%PLOTCENTRALBODY Plots a spherical central body centered at the origin.
%
% INPUT
%   body : struct with field
%       .radius  - body radius [km]
%
% OUTPUT
%   h : surface handle

[xs, ys, zs] = sphere(60);

xs = body.radius * xs;
ys = body.radius * ys;
zs = body.radius * zs;

h = surf(xs, ys, zs);
shading interp
axis equal
end