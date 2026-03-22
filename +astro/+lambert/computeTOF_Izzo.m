function [T, y, psi] = computeTOF_Izzo(x, lambda, M)
%COMPUTETOF_IZZO Non-dimensional time-of-flight for Izzo/Lancaster-Blanchard form.
%
% INPUTS
%   x      : Lancaster-Blanchard variable
%   lambda : geometry parameter in [-1,1]
%   M      : number of revolutions (here use M = 0)
%
% OUTPUTS
%   T      : non-dimensional time of flight
%   y      : auxiliary variable
%   psi    : auxiliary angle/anomaly
%
% Notes:
%   Valid for elliptic/parabolic/hyperbolic cases in the 0-rev setting.
%   This is a direct implementation of the x-lambda TOF structure used by Izzo.

if nargin < 3
    M = 0;
end

y = sqrt(1 - lambda^2 * (1 - x^2));

if abs(1 - x) < 1e-8
    % Parabolic limit from Izzo paper
    T = (2/3) * (1 - lambda^3);
    psi = 0;
    return
end

if x < 1
    % Elliptic case
    arg = x*y + lambda*(1 - x^2);
    arg = max(-1, min(1, arg));
    psi = acos(arg);

    % Correct single-rev TOF formula
    T = ((psi + M*pi) / sqrt(1 - x^2) - x + lambda*y) / (1 - x^2);
else
    % Hyperbolic case
    arg = x*y - lambda*(x^2 - 1);
    if arg < 1
        arg = 1;
    end
    psi = acosh(arg);

    T = ((lambda*y - x) - psi / sqrt(x^2 - 1)) / (x^2 - 1);
end
end