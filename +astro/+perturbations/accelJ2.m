function aJ2 = accelJ2(r, mu, J2, Re)
%ACCELJ2 J2 perturbation acceleration.
%
% INPUTS
%   r  : position vector [3x1 km]
%   mu : gravitational parameter [km^3/s^2]
%   J2 : second zonal harmonic [-]
%   Re : equatorial radius [km]
%
% OUTPUT
%   aJ2 : J2 acceleration [3x1 km/s^2]

r = r(:);
x = r(1); y = r(2); z = r(3);

rmag = norm(r);
if rmag <= 0
    error('Position norm must be positive.');
end

z2 = z^2;
r2 = rmag^2;

fac = 1.5 * J2 * mu * Re^2 / rmag^5;

aJ2 = fac * [x * (5*z2/r2 - 1);
             y * (5*z2/r2 - 1);
             z * (5*z2/r2 - 3)];
end