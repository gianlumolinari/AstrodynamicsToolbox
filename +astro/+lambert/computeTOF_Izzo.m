function out = computeTOF_Izzo(x, lambda, M)
%COMPUTETOF_IZZO Izzo/Lancaster-Blanchard non-dimensional TOF and derivatives.
%
% INPUTS
%   x      : Lancaster-Blanchard variable
%   lambda : geometry parameter in [-1,1]
%   M      : number of revolutions (use M = 0 here)
%
% OUTPUT
%   out : struct with fields
%       .T
%       .y
%       .psi
%       .dTdx
%       .d2Tdx2
%       .d3Tdx3
%
% Notes:
%   Based on Izzo (2015), Eq. (18), Eq. (21), Eq. (22), Eq. (23).
%   Here we use the 0-revolution case, but the derivative formulas are
%   written with M retained in the TOF expression.

if nargin < 3
    M = 0;
end

y = sqrt(1 - lambda^2 * (1 - x^2));

% Parabolic limit
if abs(x - 1) < 1e-10
    T = (2/3) * (1 - lambda^3);
    dTdx = (2/5) * (lambda^5 - 1);
    
    % Higher derivatives are not used exactly at x = 1 in practice.
    d2Tdx2 = NaN;
    d3Tdx3 = NaN;
    
    out.T = T;
    out.y = y;
    out.psi = 0;
    out.dTdx = dTdx;
    out.d2Tdx2 = d2Tdx2;
    out.d3Tdx3 = d3Tdx3;
    return
end

if x < 1
    % Elliptic case
    arg = x*y + lambda*(1 - x^2);
    arg = max(-1, min(1, arg));
    psi = acos(arg);
    
    T = ((psi + M*pi) / sqrt(1 - x^2) - x + lambda*y) / (1 - x^2);
else
    % Hyperbolic case
    arg = x*y - lambda*(x^2 - 1);
    arg = max(1, arg);
    psi = acosh(arg);
    
    T = ((lambda*y - x) - psi / sqrt(x^2 - 1)) / (x^2 - 1);
end

% Analytical derivatives from Izzo paper, Eq. (22)
omx2 = 1 - x^2;

% Guard against singular evaluations extremely close to x = ±1
if abs(omx2) < 1e-12 || abs(y) < 1e-12
    % fallback numerical derivatives
    h = 1e-6;
    Tp = computeTOF_Izzo(x + h, lambda, M);
    Tm = computeTOF_Izzo(x - h, lambda, M);
    dTdx = (Tp.T - Tm.T) / (2*h);
    d2Tdx2 = (Tp.T - 2*T + Tm.T) / (h^2);
    
    Tpp = computeTOF_Izzo(x + 2*h, lambda, M);
    Tmm = computeTOF_Izzo(x - 2*h, lambda, M);
    d3Tdx3 = (Tpp.T - 2*Tp.T + 2*Tm.T - Tmm.T) / (2*h^3);
else
    dTdx = (3*T*x - 2 + 2*lambda^3 * x / y) / omx2;
    d2Tdx2 = (3*T + 5*x*dTdx + 2*(1 - lambda^2)*lambda^3 / y^3) / omx2;
    d3Tdx3 = (7*x*d2Tdx2 + 8*dTdx - 6*(1 - lambda^2)*lambda^5 * x / y^5) / omx2;
end

out.T = T;
out.y = y;
out.psi = psi;
out.dTdx = dTdx;
out.d2Tdx2 = d2Tdx2;
out.d3Tdx3 = d3Tdx3;
end