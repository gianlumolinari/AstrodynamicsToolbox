function plotFlyby2D(vInfIn, vInfOutReq, mu, rp, turnSign, bodyRadius, bodyName)
%PLOTFLYBY2D Simple 2D body-centered flyby geometry plot.
%
% INPUTS
%   vInfIn      : incoming v_infinity vector [3x1 km/s]
%   vInfOutReq  : required outgoing v_infinity vector [3x1 km/s]
%   mu          : body gravitational parameter [km^3/s^2]
%   rp          : chosen periapsis radius [km]
%   turnSign    : +1 or -1
%   bodyRadius  : body radius [km]
%   bodyName    : name for title

if nargin < 7
    bodyName = 'Flyby Body';
end

fly = astro.mga.flybyTurn(vInfIn, mu, rp, turnSign);
vInfOutAch = fly.vInfOut;

scale = 15000 / max([norm(vInfIn), norm(vInfOutReq), norm(vInfOutAch)]);

vin  = vInfIn(:) * scale;
vreq = vInfOutReq(:) * scale;
vach = vInfOutAch(:) * scale;

figure('Color','w');
hold on
axis equal
grid on

% Body disk
ang = linspace(0, 2*pi, 300);
fill(bodyRadius*cos(ang), bodyRadius*sin(ang), [0.85 0.92 1.0], ...
    'EdgeColor', 'k', 'LineWidth', 1.2);

% Periapsis circle
plot(rp*cos(ang), rp*sin(ang), '--', 'LineWidth', 1.2);

% Velocity vectors
quiver(0,0,vin(1),vin(2),0,'LineWidth',2,'MaxHeadSize',0.5);
quiver(0,0,vreq(1),vreq(2),0,'LineWidth',2,'MaxHeadSize',0.5);
quiver(0,0,vach(1),vach(2),0,'LineWidth',2,'MaxHeadSize',0.5);

xlabel('x [km]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('y [km]', 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('%s Flyby Geometry', upper(bodyName)), ...
    'FontSize', 15, 'FontWeight', 'bold');

legend(bodyName, 'Periapsis radius', 'Incoming v_\infty^-', ...
    'Required outgoing v_\infty^+', 'Achievable outgoing v_\infty^+', ...
    'Location', 'best');

text(0.05*rp, -1.35*rp, sprintf('Turn angle = %.3f deg', fly.deltaDeg), ...
    'FontSize', 11, 'FontWeight', 'bold');