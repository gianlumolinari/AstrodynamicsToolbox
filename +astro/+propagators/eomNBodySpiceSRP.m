function dx = eomNBodySpiceSRP(t, x, config)
%EOMNBODYSPICESRP High-fidelity heliocentric EOM using SPICE ephemerides.
%
% INPUT
%   t      : time since initial epoch [s]
%   x      : state [rx; ry; rz; vx; vy; vz] wrt Sun in J2000 [km, km/s]
%   config : struct with fields:
%       .et0         initial SPICE ET [s]
%       .muSun       Sun GM [km^3/s^2]
%       .bodyNames   cell array of perturbing bodies
%       .muBodies    corresponding GMs [km^3/s^2]
%       .useSRP      true/false
%       .Cr          reflectivity coefficient
%       .A_over_m    area-to-mass ratio [m^2/kg]
%
% OUTPUT
%   dx     : state derivative

    r = x(1:3);
    v = x(4:6);

    rnorm = norm(r);
    a = -config.muSun * r / rnorm^3;

    et = config.et0 + t;

    for k = 1:numel(config.bodyNames)
        stBody = cspice_spkezr(config.bodyNames{k}, et, 'J2000', 'NONE', 'SUN');
        rBody = stBody(1:3);

        rel = rBody - r;
        a = a + config.muBodies(k) * ( rel / norm(rel)^3 - rBody / norm(rBody)^3 );
    end

    if isfield(config,'useSRP') && config.useSRP
        a = a + astro.propagators.accelSRP(r, config.Cr, config.A_over_m);
    end

    dx = [v; a];
end