function h = plotCentralBody(body, faceColor)
%PLOTCENTRALBODY Plots a spherical central body centered at the origin.
%
% INPUT
%   body : struct with field
%       .radius  - body radius [km]
%   faceColor : optional RGB color, default soft Earth-like blue
%
% OUTPUT
%   h : surface handle

if nargin < 2 || isempty(faceColor)
    faceColor = [0.78 0.86 0.96];
end

[xs, ys, zs] = sphere(80);

xs = body.radius * xs;
ys = body.radius * ys;
zs = body.radius * zs;

h = surf(xs, ys, zs, ...
    'FaceColor', faceColor, ...
    'EdgeColor', 'none', ...
    'FaceLighting', 'gouraud', ...
    'AmbientStrength', 0.35, ...
    'DiffuseStrength', 0.75, ...
    'SpecularStrength', 0.15);

axis equal
shading interp
camlight headlight
material dull
end