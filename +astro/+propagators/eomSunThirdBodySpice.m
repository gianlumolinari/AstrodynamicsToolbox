function dx = eomSunThirdBodySpice(t, x, et0, bodyNames, muBodies)
%EOMSUNTHIRDBODYSPICE Sun-centered EOM with SPICE third-body perturbations.
%
% INPUTS
%   t         : integration time since initial epoch [s]
%   x         : state vector [rx; ry; rz; vx; vy; vz] in heliocentric J2000
%   et0       : initial SPICE ephemeris time [s past J2000]
%   bodyNames : cell array of perturbing body names, e.g. {'EARTH','MARS','JUPITER BARYCENTER'}
%   muBodies  : vector of gravitational parameters [km^3/s^2]
%
% OUTPUT
%   dx        : state derivative [vx; vy; vz; ax; ay; az]
%
% MODEL
%   Heliocentric equations of motion with:
%     - Sun as central body
%     - third-body perturbations using SPICE body positions wrt SUN
%
%   a = -muSun*r/|r|^3 + sum_i mu_i * [ (r_i-r)/|r_i-r|^3 - r_i/|r_i|^3 ]

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

    dx = [v; a];
end