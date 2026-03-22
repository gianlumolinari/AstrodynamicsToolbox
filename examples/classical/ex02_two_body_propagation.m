clc;
clear;
close all;

startup

% ==========================================================
% Example: two-body orbit propagation around Earth
% ==========================================================

earth = astro.bodies.earth();

% Classical orbital elements
a     = 7000;              % km
e     = 0.01;              % -
i     = deg2rad(30);       % rad
RAAN  = deg2rad(40);       % rad
omega = deg2rad(60);       % rad
theta = deg2rad(20);       % rad

% Convert to Cartesian state
[r0, v0] = astro.coords.coe2rv(a, e, i, RAAN, omega, theta, earth.mu);
x0 = [r0; v0];

% One orbital period
T = 2*pi*sqrt(a^3/earth.mu);
tspan = [0 T];

% Propagation options
opts.RelTol = 1e-12;
opts.AbsTol = 1e-12;
opts.Solver = 'ode113';

% Propagate
out = astro.propagators.propagate( ...
    @(t,x) astro.propagators.eomTwoBody(t, x, earth.mu), ...
    tspan, x0, opts);

% Compute diagnostics
N = size(out.x,1);
epsHist = zeros(N,1);
hHist   = zeros(N,3);

for k = 1:N
    r = out.x(k,1:3).';
    v = out.x(k,4:6).';
    epsHist(k) = astro.utils.specificEnergy(r, v, earth.mu);
    hHist(k,:) = astro.utils.angularMomentum(r, v).';
end

energyDrift = max(abs(epsHist - epsHist(1)));
hDrift      = max(vecnorm(hHist - hHist(1,:), 2, 2));

fprintf('Max energy drift            = %.6e km^2/s^2\n', energyDrift);
fprintf('Max angular momentum drift  = %.6e km^2/s\n', hDrift);

% Plot
figure
astro.plot.plotCentralBody(earth);
hold on
astro.plot.plotOrbit3D(out.x(:,1:3), 'LineWidth', 1.5);
title('Two-body orbit propagation around Earth')
legend('Earth','Orbit','Location','best')