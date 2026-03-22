function body = mars()
%MARS Returns basic physical data for Mars.
%
% OUTPUT
%   body : struct with fields
%       .name   - body name
%       .mu     - gravitational parameter [km^3/s^2]
%       .radius - mean equatorial radius [km]
%       .J2     - second zonal harmonic [-]

body.name   = 'Mars';
body.mu     = 42828.375214;   % km^3/s^2
body.radius = 3396.19;        % km
body.J2     = 1.96045e-3;     % -
end
