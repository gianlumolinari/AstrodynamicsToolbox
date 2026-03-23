function out = isInShadowConical(rSc, rSun, rBody, bodyRadius, sunRadius)
%ISINSHADOWCONICAL Conical umbra/penumbra eclipse test.
%
% INPUTS
%   rSc        : spacecraft position [3x1 km]
%   rSun       : Sun position [3x1 km]
%   rBody      : occulting body center [3x1 km]
%   bodyRadius : occulting body radius [km]
%   sunRadius  : Sun radius [km]
%
% OUTPUT
%   out : struct with fields
%       .inUmbra
%       .inPenumbra
%       .inShadow
%       .state          'sunlight' | 'penumbra' | 'umbra'
%       .alongTrack
%       .crossTrack
%       .rUmbra
%       .rPenumbra
%
% NOTES
%   Uses straight-line conical geometry with the Sun treated as a sphere.

rSc = rSc(:);
rSun = rSun(:);
rBody = rBody(:);

s = rSun - rBody;
D = norm(s);

if D <= bodyRadius || D <= sunRadius
    error('Invalid body/Sun geometry for eclipse cones.');
end

sHat = s / D;
aHat = -sHat;  % shadow axis

rho = rSc - rBody;
alongTrack = dot(rho, aHat);
rhoPerp = rho - alongTrack * aHat;
crossTrack = norm(rhoPerp);

% Umbra and penumbra cone lengths
Lumbra = D * bodyRadius / (sunRadius - bodyRadius);
Lpen   = D * bodyRadius / (sunRadius + bodyRadius);

% Default
rUmbra = NaN;
rPenumbra = NaN;
inUmbra = false;
inPenumbra = false;

if alongTrack > 0
    % Cone radii at the spacecraft axial distance
    rUmbra = bodyRadius * (1 - alongTrack / Lumbra);
    rPenumbra = bodyRadius * (1 + alongTrack / Lpen);

    if alongTrack < Lumbra && rUmbra > 0 && crossTrack <= rUmbra
        inUmbra = true;
    elseif crossTrack <= rPenumbra
        inPenumbra = true;
    end
end

inShadow = inUmbra || inPenumbra;

if inUmbra
    state = 'umbra';
elseif inPenumbra
    state = 'penumbra';
else
    state = 'sunlight';
end

out.inUmbra = inUmbra;
out.inPenumbra = inPenumbra;
out.inShadow = inShadow;
out.state = state;
out.alongTrack = alongTrack;
out.crossTrack = crossTrack;
out.rUmbra = rUmbra;
out.rPenumbra = rPenumbra;
out.Lumbra = Lumbra;
out.Lpenumbra = Lpen;
end