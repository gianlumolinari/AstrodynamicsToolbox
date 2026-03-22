function coe = rv2coe(r, v, mu)
%RV2COE Converts Cartesian state vectors to classical orbital elements.
%
% INPUTS
%   r   : position vector [km]
%   v   : velocity vector [km/s]
%   mu  : gravitational parameter [km^3/s^2]
%
% OUTPUT
%   coe : struct with fields
%       .a
%       .e
%       .i
%       .RAAN
%       .omega
%       .theta

r = r(:);
v = v(:);

R = norm(r);
V = norm(v);

h_vec = cross(r, v);
h = norm(h_vec);

k_vec = [0; 0; 1];
n_vec = cross(k_vec, h_vec);
n = norm(n_vec);

e_vec = (1/mu) * ((V^2 - mu/R)*r - dot(r,v)*v);
e = norm(e_vec);

energy = V^2/2 - mu/R;

if abs(e - 1) > 1e-12
    a = -mu / (2*energy);
else
    a = inf;
end

i = acos(h_vec(3)/h);

% RAAN
if n > 1e-12
    RAAN = acos(n_vec(1)/n);
    if n_vec(2) < 0
        RAAN = 2*pi - RAAN;
    end
else
    RAAN = 0;
end

% Argument of periapsis
if n > 1e-12 && e > 1e-12
    omega = acos(dot(n_vec,e_vec)/(n*e));
    if e_vec(3) < 0
        omega = 2*pi - omega;
    end
else
    omega = 0;
end

% True anomaly
if e > 1e-12
    theta = acos(dot(e_vec,r)/(e*R));
    if dot(r,v) < 0
        theta = 2*pi - theta;
    end
else
    theta = 0;
end

coe.a     = a;
coe.e     = e;
coe.i     = i;
coe.RAAN  = RAAN;
coe.omega = omega;
coe.theta = theta;
end