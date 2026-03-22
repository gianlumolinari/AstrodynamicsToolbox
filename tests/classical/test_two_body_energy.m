clc;
clear;
close all;

startup

earth = astro.bodies.earth();

% Test orbit
a     = 7000;              
e     = 0.01;
i     = deg2rad(30);
RAAN  = deg2rad(40);
omega = deg2rad(60);
theta = deg2rad(20);

[r0, v0] = astro.coords.coe2rv(a, e, i, RAAN, omega, theta, earth.mu);
x0 = [r0; v0];

T = 2*pi*sqrt(a^3/earth.mu);
tspan = [0 T];

opts.RelTol = 1e-12;
opts.AbsTol = 1e-12;
opts.Solver = 'ode113';

out = astro.propagators.propagate( ...
    @(t,x) astro.propagators.eomTwoBody(t, x, earth.mu), ...
    tspan, x0, opts);

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

energyTol = 1e-9;
hTol      = 1e-6;

if energyDrift < energyTol && hDrift < hTol
    disp('PASS: Two-body propagation conserves energy and angular momentum.')
else
    error('FAIL: Two-body propagation did not satisfy conservation tolerances.')
end