function aDrag = accelDrag(r, v, Cd, AoverM, rho0, H, Re, omegaBody)
%ACCELDRAG Simple exponential drag model.
%
% INPUTS
%   r         : position vector [3x1 km]
%   v         : inertial velocity vector [3x1 km/s]
%   Cd        : drag coefficient [-]
%   AoverM    : area-to-mass ratio [m^2/kg]
%   rho0      : reference density [kg/m^3]
%   H         : scale height [km]
%   Re        : body radius [km]
%   omegaBody : body rotation rate [rad/s]
%
% OUTPUT
%   aDrag : drag acceleration [3x1 km/s^2]
%
% NOTES
%   Internally converts to SI for drag computation, then returns km/s^2.

r = r(:);
v = v(:);

rmag = norm(r);
h = rmag - Re;

if h < 0
    error('Spacecraft is below the reference radius.');
end

rho = rho0 * exp(-h / H);

omegaVec = [0; 0; omegaBody];
vAtm = cross(omegaVec, r);      % km/s
vRel = v - vAtm;                % km/s
vRel_m = 1000 * vRel;           % m/s
vRelMag = norm(vRel_m);

if vRelMag == 0
    aDrag = zeros(3,1);
    return
end

aDrag_m = -0.5 * Cd * AoverM * rho * vRelMag * vRel_m;  % m/s^2
aDrag = aDrag_m / 1000;                                 % km/s^2
end