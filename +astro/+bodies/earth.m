function body = earth()
%EARTH Returns basic physical data for Earth.
%
% OUTPUT
%   body : struct with fields
%       .name   - body name
%       .mu     - gravitational parameter [km^3/s^2]
%       .radius - mean equatorial radius [km]
%       .J2     - second zonal harmonic [-]

body.name   = 'Earth';
body.mu     = 398600.4418;      % km^3/s^2
body.radius = 6378.1363;        % km
body.J2     = 1.08262668e-3;    % -
end
