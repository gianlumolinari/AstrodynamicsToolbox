function out = gravityAssist3D(bodyName, rp, thetaB, vInfIn, vPlanet, et)
%GRAVITYASSIST3D 3D patched-conic gravity-assist map in J2000.
%
% INPUTS
%   bodyName : string, e.g. 'earth', 'venus', 'jupiter'
%   rp       : periapsis radius from body center [km]
%   thetaB   : B-plane angle [rad]
%   vInfIn   : incoming hyperbolic excess velocity wrt body [3x1 km/s]
%   vPlanet  : heliocentric body velocity in J2000 [3x1 km/s]
%   et       : SPICE ephemeris time [s past J2000]
%
% OUTPUT
%   out : struct with fields
%       .bodyName
%       .mu
%       .rp
%       .thetaBa
%       .vInfIn
%       .vInfOut
%       .vInfMagIn
%       .vInfMagOut
%       .deltaRad
%       .deltaDeg
%       .Bmag
%       .Bhat
%       .Shat
%       .That
%       .Rhat
%       .vOutHelio
%
% NOTES
%   - Uses the planet pole from IAU_<BODY> to define the B-plane frame.
%   - The outgoing heliocentric velocity is
%         vOutHelio = vPlanet + vInfOut
%   - This is a patched-conic flyby map. It does not model the passage
%     through the SOI with full n-body dynamics.

    body = astro.bodies.getBody(lower(bodyName));
    mu = body.mu;

    vInfIn = vInfIn(:);
    vPlanet = vPlanet(:);

    vInfMag = norm(vInfIn);

    if vInfMag <= 0
        error('vInfIn must have nonzero magnitude.');
    end
    if rp <= 0
        error('rp must be positive.');
    end

    % Planet pole in J2000 from body-fixed frame
    frameBody = ['IAU_' upper(bodyName)];
    M = cspice_pxform(frameBody, 'J2000', et);
    K = M * [0; 0; 1];
    K = K / norm(K);

    % Incoming asymptote direction
    S = vInfIn / vInfMag;

    % Build B-plane frame
    c = cross(S, K);
    if norm(c) < 1e-12
        error('Incoming v_inf is nearly parallel to the planet pole. B-plane ill-defined.');
    end

    T = c / norm(c);
    R = cross(S, T);
    R = R / norm(R);

    % Hyperbolic turning angle
    ecc = 1 + rp * vInfMag^2 / mu;
    delta = 2 * asin(1 / ecc);

    % Equivalent B magnitude
    Bmag = sqrt(rp^2 + 2*mu*rp / vInfMag^2);

    % B-plane direction from thetaB
    Bhat = cos(thetaB) * T + sin(thetaB) * R;
    Bhat = Bhat / norm(Bhat);

    % Rodrigues rotation of incoming v_inf about Bhat
    vInfOut = vInfIn * cos(delta) ...
            + cross(Bhat, vInfIn) * sin(delta) ...
            + Bhat * (dot(Bhat, vInfIn)) * (1 - cos(delta));

    % Heliocentric outgoing velocity
    vOutHelio = vPlanet + vInfOut;

    out.bodyName   = lower(bodyName);
    out.mu         = mu;
    out.rp         = rp;
    out.thetaB     = thetaB;
    out.vInfIn     = vInfIn;
    out.vInfOut    = vInfOut;
    out.vInfMagIn  = norm(vInfIn);
    out.vInfMagOut = norm(vInfOut);
    out.deltaRad   = delta;
    out.deltaDeg   = rad2deg(delta);
    out.Bmag       = Bmag;
    out.Bhat       = Bhat;
    out.Shat       = S;
    out.That       = T;
    out.Rhat       = R;
    out.vOutHelio  = vOutHelio;
end