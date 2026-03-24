function out = propagateToSection(x0, tMax, mu, xSection, direction, stopOnFirst)
%PROPAGATETOSECTION Propagate a CR3BP trajectory until it crosses x = xSection.
%
% INPUTS
%   x0          : initial state [6x1]
%   tMax        : maximum propagation time
%                 use positive for forward, negative for backward
%   mu          : CR3BP mass parameter
%   xSection    : section plane x = constant
%   direction   : event crossing direction
%                 +1 = increasing through section
%                 -1 = decreasing through section
%                  0 = any crossing
%   stopOnFirst : true/false, stop at first crossing
%
% OUTPUT
%   out : struct with fields
%       .t
%       .x
%       .tCross
%       .xCross
%       .hit

if nargin < 6 || isempty(stopOnFirst)
    stopOnFirst = true;
end
if nargin < 5 || isempty(direction)
    direction = 0;
end

x0 = x0(:);

opts = odeset( ...
    'RelTol', 1e-12, ...
    'AbsTol', 1e-12, ...
    'Events', @(t,x) localSectionEvent(t, x, xSection, direction, stopOnFirst));

[t, x, tE, xE] = ode113( ...
    @(t,x) astro.cr3bp.eomCR3BP(t, x, mu), ...
    [0 tMax], x0, opts);

out.t = t;
out.x = x;
out.hit = ~isempty(tE);

if out.hit
    out.tCross = tE(1);
    out.xCross = xE(1,:).';
else
    out.tCross = NaN;
    out.xCross = NaN(6,1);
end

end

function [value, isterminal, dir] = localSectionEvent(t, x, xSection, direction, stopOnFirst)
% avoid immediate trigger at t = 0 if starting very near section
if abs(t) < 1e-12
    value = 1;
else
    value = x(1) - xSection;
end

if stopOnFirst
    isterminal = 1;
else
    isterminal = 0;
end

dir = direction;
end