function man = hohmann(mu, r1, r2)
%HOHMANN Computes the parameters of a Hohmann transfer between two circular orbits.
%
% INPUTS
%   mu : gravitational parameter [km^3/s^2]
%   r1 : initial circular orbit radius [km]
%   r2 : final circular orbit radius [km]
%
% OUTPUT
%   man : struct with fields
%       .r1
%       .r2
%       .aTrans
%       .vCirc1
%       .vCirc2
%       .vPeriTrans
%       .vApoTrans
%       .dv1
%       .dv2
%       .dvTot
%       .tof

aTrans = 0.5 * (r1 + r2);

vCirc1 = sqrt(mu / r1);
vCirc2 = sqrt(mu / r2);

vPeriTrans = sqrt(mu * (2/r1 - 1/aTrans));
vApoTrans  = sqrt(mu * (2/r2 - 1/aTrans));

dv1 = abs(vPeriTrans - vCirc1);
dv2 = abs(vCirc2 - vApoTrans);
dvTot = dv1 + dv2;

tof = pi * sqrt(aTrans^3 / mu);

man.r1 = r1;
man.r2 = r2;
man.aTrans = aTrans;
man.vCirc1 = vCirc1;
man.vCirc2 = vCirc2;
man.vPeriTrans = vPeriTrans;
man.vApoTrans = vApoTrans;
man.dv1 = dv1;
man.dv2 = dv2;
man.dvTot = dvTot;
man.tof = tof;
end