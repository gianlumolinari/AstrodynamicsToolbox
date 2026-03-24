function aSRP = accelSRP(rSunToSc, Cr, A_over_m)
%ACCELSRP Solar radiation pressure acceleration in heliocentric frame.
%
% INPUTS
%   rSunToSc : spacecraft position wrt Sun [km]
%   Cr       : reflectivity coefficient [-]
%   A_over_m : area-to-mass ratio [m^2/kg]
%
% OUTPUT
%   aSRP     : SRP acceleration [km/s^2]
%
% MODEL
%   a = P0 * Cr * (A/m) * (AU/r)^2 * rhat
%
% where:
%   P0  = solar radiation pressure at 1 AU [N/m^2]
%   AU  = astronomical unit [m]
%   r   = Sun-spacecraft distance [m]
%
% Direction is radially away from the Sun.

    P0 = 4.56e-6;          % N/m^2 at 1 AU
    AU = 149597870700;     % m

    r_m = rSunToSc(:) * 1e3;
    rnorm_m = norm(r_m);

    if rnorm_m == 0
        error('SRP undefined at zero Sun-spacecraft distance.');
    end

    rhat = r_m / rnorm_m;

    a_mps2 = P0 * Cr * A_over_m * (AU / rnorm_m)^2 * rhat;
    aSRP = a_mps2 / 1e3;   % convert to km/s^2
end