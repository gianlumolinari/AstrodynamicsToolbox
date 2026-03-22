function h = plotOrbit3D(r, varargin)
%PLOTORBIT3D Plots a 3D trajectory.
%
% INPUT
%   r : Nx3 matrix of position history [km]
%
% OPTIONAL INPUT
%   Additional plotting options, e.g.:
%       plotOrbit3D(r, 'LineWidth', 1.5)
%
% OUTPUT
%   h : line handle

h = plot3(r(:,1), r(:,2), r(:,3), varargin{:});
grid on
axis equal
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
end