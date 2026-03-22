function [r_eci, v_eci] = coe2rv(a, e, i, RAAN, omega, theta, mu)
%COE2RV Converts classical orbital elements to Cartesian state vectors.
%
% INPUTS
%   a      : semi-major axis [km]
%   e      : eccentricity [-]
%   i      : inclination [rad]
%   RAAN   : right ascension of ascending node [rad]
%   omega  : argument of periapsis [rad]
%   theta  : true anomaly [rad]
%   mu     : gravitational parameter [km^3/s^2]
%
% OUTPUTS
%   r_eci  : position vector in inertial frame [km]
%   v_eci  : velocity vector in inertial frame [km/s]

% Semi-latus rectum
p = a * (1 - e^2);

% Position and velocity in perifocal frame
r_pf = (p / (1 + e*cos(theta))) * [cos(theta); sin(theta); 0];
v_pf = sqrt(mu / p) * [-sin(theta); e + cos(theta); 0];

% Rotation matrices
R3_W = [ cos(RAAN), -sin(RAAN), 0;
         sin(RAAN),  cos(RAAN), 0;
                 0,          0, 1];

R1_i = [1,       0,        0;
        0,  cos(i), -sin(i);
        0,  sin(i),  cos(i)];

R3_w = [ cos(omega), -sin(omega), 0;
         sin(omega),  cos(omega), 0;
                  0,           0, 1];

% Perifocal to inertial rotation
Q_pX = R3_W * R1_i * R3_w;

% Transform to inertial frame
r_eci = Q_pX * r_pf;
v_eci = Q_pX * v_pf;
end