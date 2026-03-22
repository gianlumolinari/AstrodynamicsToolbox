function body = moon()
%MOON Returns basic physical data for the Moon.
%
% OUTPUT
%   body : struct with fields
%       .name   - body name
%       .mu     - gravitational parameter [km^3/s^2]
%       .radius - mean radius [km]
%       .J2     - second zonal harmonic [-]

body.name   = 'Moon';
body.mu     = 4902.800066;   % km^3/s^2
body.radius = 1737.4;        % km
body.J2     = 0.0;           % -
end
