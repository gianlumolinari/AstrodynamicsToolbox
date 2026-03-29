function mu = getSystemMu(systemName)
%GETSYSTEMMU Return CR3BP mass parameter for supported systems.

switch lower(string(systemName))
    case "earth_moon"
        mu = 0.012150585609624;

    case "sun_earth"
        mu = 3.040357143198e-6;

    otherwise
        error('Unsupported CR3BP system: %s', systemName);
end
end