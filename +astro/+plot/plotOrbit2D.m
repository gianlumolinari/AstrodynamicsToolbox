function h = plotOrbit2D(r, varargin)
%PLOTORBIT2D Plots a 2D trajectory.
%
% INPUT
%   r : Nx2 or Nx3 matrix of position history [km]
%
% OPTIONAL INPUT
%   Additional plotting options
%
% OUTPUT
%   h : line handle

if size(r,2) < 2
    error('plotOrbit2D requires an Nx2 or Nx3 position array.');
end

h = plot(r(:,1), r(:,2), varargin{:});
grid on
axis equal
xlabel('x [km]')
ylabel('y [km]')
end