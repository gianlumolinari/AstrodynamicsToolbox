function a3b = accelThirdBody(rSc, rThirdBody, muThird)
%ACCELTHIRDBODY Third-body perturbation acceleration.
%
% INPUTS
%   rSc        : spacecraft position relative to central body [3x1 km]
%   rThirdBody : third-body position relative to central body [3x1 km]
%   muThird    : third-body gravitational parameter [km^3/s^2]
%
% OUTPUT
%   a3b        : third-body perturbation acceleration [3x1 km/s^2]
%
% FORMULA
%   a3b = mu3 * ( (r3 - r)/|r3-r|^3 - r3/|r3|^3 )

rSc = rSc(:);
rThirdBody = rThirdBody(:);

rho = rThirdBody - rSc;

rhoMag = norm(rho);
r3Mag = norm(rThirdBody);

if rhoMag <= 0 || r3Mag <= 0
    error('Invalid geometry in accelThirdBody.');
end

a3b = muThird * (rho / rhoMag^3 - rThirdBody / r3Mag^3);
end