function h = plotTransfer2D(r1, r2, varargin)
%PLOTTRANSFER2D Plots a straight-line geometry between two position vectors.
%
% INPUTS
%   r1 : initial position vector [2x1 or 3x1]
%   r2 : final position vector [2x1 or 3x1]
%
% OPTIONAL INPUT
%   Additional plotting options
%
% OUTPUT
%   h : line handle
%
% Note:
%   This utility is for simple geometry visualization. It does not propagate
%   the Lambert arc. It just draws the segment between r1 and r2.

r1 = r1(:);
r2 = r2(:);

h = plot([r1(1) r2(1)], [r1(2) r2(2)], varargin{:});
grid on
axis equal
xlabel('x [km]')
ylabel('y [km]')
end