function dx = eomSunThirdBodySpiceSRP(t, x, et0, bodyNames, muBodies, Cr, A_over_m)
%EOMSUNTHIRDBODYSPICESRP Sun-centered EOM with SPICE third bodies + SRP.
%
% INPUTS
%   t         : integration time since initial epoch [s]
%   x         : state [rx; ry; rz; vx; vy; vz] wrt Sun in J2000 [km, km/s]
%   et0       : initial SPICE ephemeris time [s past J2000]
%   bodyNames : cell array of perturbing body names
%   muBodies  : corresponding gravitational parameters [km^3/s^2]
%   Cr        : SRP reflectivity coefficient [-]
%   A_over_m  : area-to-mass ratio [m^2/kg]
%
% OUTPUT
%   dx        : state derivative

    r = x(1:3);
    v = x(4:6);

    sun = astro.bodies.getBody('sun');
    muSun = sun.mu;

    rnorm = norm(r);
    a = -muSun * r / rnorm^3;

    et = et0 + t;

    for k = 1:numel(bodyNames)
        stateBody = cspice_spkezr(bodyNames{k}, et, 'J2000', 'NONE', 'SUN');
        rBody = stateBody(1:3);

        rel = rBody - r;
        a = a + muBodies(k) * ( rel / norm(rel)^3 - rBody / norm(rBody)^3 );
    end

    a = a + astro.propagators.accelSRP(r, Cr, A_over_m);

    dx = [v; a];
end