function body = sun()
%SUN Returns basic physical data for the Sun.
%
% OUTPUT
%   body : struct with fields
%       .name   - body name
%       .mu     - gravitational parameter [km^3/s^2]
%       .radius - mean radius [km]
%       .J2     - second zonal harmonic [-]

body.name   = 'Sun';
body.mu     = 1.32712440018e11;   % km^3/s^2
body.radius = 695700;             % km
body.J2     = 0.0;                % -
end
